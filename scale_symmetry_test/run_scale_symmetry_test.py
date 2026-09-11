"""cheaseBS scale-mapping and response-symmetry test — one sweep, two questions.

Runs on NERSC (needs the compiled CHEASE binary via TPED/config/user_config.yaml).

Both open questions from the 2026-08-18 five-point 132543 run are properties of
the same curve, so one symmetric sweep answers both.

QUESTION 1 — scale mapping.
    `apply_mtanh_ped(scale_height=s)` multiplies the mtanh STEP amplitude a_y0,
    not the pedestal-top value, which is y_sep + a_y0 + core. Measured on the
    132543 seed, s=1.3 gave +14.4% at the pedestal top for ne and +21.3% for Te
    — so a scan box quoted in scale factors is much narrower than the same
    numbers read as pedestal-top fractions, and the survey's absolute bounds
    table (T_e,ped 0.2-0.8 keV, n_e,ped 3-7e19) cannot be mapped onto scale
    factors without this curve. Is the mapping linear? What is the slope per
    variable?

QUESTION 2 — response symmetry.
    ne x1.3 took 19 cheaseBS iterations and moved Ip error 1.62% -> 0.77%.
    ne x0.7 took 2 iterations and moved it 1.62% -> 1.62%, essentially not at
    all. The profile input is exactly symmetric (+/-14.4% at the pedestal top),
    so the asymmetry is in cheaseBS's response. Candidates: the floored p_fast
    split, quasineutrality clipping nz at low density, or qspec pinning the
    current asymmetrically. Which?

The sweep is symmetric in scale factor about 1.0 and covers both variables, so
each row carries the measured pedestal-top ratio (Q1) alongside the iteration
count, Ip error and bootstrap profile (Q2).

Usage (NERSC):
    python run_scale_symmetry_test.py --dirpath /path/to/132543
    python run_scale_symmetry_test.py --dirpath ... --scales 0.7,0.85,1.0,1.15,1.3
    python run_scale_symmetry_test.py --dirpath ... --vars ne        # ne only
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime

import numpy as np

from TPED.projects.discharge_tools.src.discharge_data import DischargeData
from TPED.projects.discharge_tools.src.discharge_physics import DischargePhysics
from TPED.projects.discharge_tools.src.cheasebs_runner import (
    CheasebsAcceptance, read_acceptance)
from TPED.projects.discharge_tools.src.transforms.mtanh_transforms import fit_mtanh

DEFAULT_DIRPATH = "/global/homes/j/joeschm/data/ST_research/NSTXU_discharges/132543"

# Symmetric about 1.0 in scale factor: the whole point is comparing +s against
# -s at equal magnitude, so any asymmetry in the result is cheaseBS's.
DEFAULT_SCALES = (0.7, 0.85, 1.0, 1.15, 1.3)

# The radii the 132543 GENE scans use. q is scored here by the gate.
GENE_RADII = (0.736, 0.825)

# Matches the sparse-scan campaign so the numbers transfer directly.
CHEASEBS_OVERRIDES = {
    "max_iter": 25,
    "tol_bs": 1e-4,
    "tol_q": 1e-4,
    "tol_ip_rel": 0.02,
    "bootstrap_mix": 0.1,
    "istar_mix": 0.05,
    "plot_errors": True,
}

PED = (0.90, 0.98)          # pedestal-top window for the measured ratio
MTANH_FIT_KWARGS = dict(pedestal_weight=8.0)


def ped_top(x, y, lo=PED[0], hi=PED[1]):
    """Max of y over the pedestal-top window — the quantity bounds are quoted in."""
    m = (x >= lo) & (x <= hi)
    return float(np.max(y[m])) if m.any() else float("nan")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dirpath", default=DEFAULT_DIRPATH)
    ap.add_argument("--outroot", default=os.path.join(os.path.dirname(__file__), "runs"))
    ap.add_argument("--scales", default=",".join(str(s) for s in DEFAULT_SCALES),
                    help="comma-separated scale factors, symmetric about 1.0")
    ap.add_argument("--vars", default="ne,Te",
                    help="comma-separated variables to sweep")
    args = ap.parse_args()

    scales = [float(s) for s in args.scales.split(",")]
    variables = [v.strip() for v in args.vars.split(",")]
    # Runs nest under a date directory so an outroot accumulating many
    # campaigns stays browsable: <outroot>/<YYYYMMDD>/<YYYYMMDD_HH-MM-SS>/
    now = datetime.now()
    stamp = now.strftime("%Y%m%d_%H-%M-%S")
    workdir = os.path.abspath(
        os.path.join(args.outroot, now.strftime("%Y%m%d"), stamp))
    os.makedirs(workdir, exist_ok=True)

    print(f"[test] workdir  {workdir}")
    print(f"[test] scales   {scales}")
    print(f"[test] vars     {variables}")

    discharge = DischargeData(input_dir=args.dirpath)
    phys0 = DischargePhysics(discharge)
    x0 = phys0.rhot.values if hasattr(phys0.rhot, "values") else np.asarray(phys0.rhot)

    # Nominal fits, built once and reused for every point, so every scale factor
    # perturbs the same reference parameterization (fit noise must not enter the
    # axis definition — the same rule the campaign follows).
    fits, nominal = {}, {}
    for var in variables:
        profile, record = fit_mtanh(phys0.ds, var, **MTANH_FIT_KWARGS)
        fits[var] = profile
        nominal[var] = record["fit_params"]
        print(f"[test] {var} fit rms_relative = {record['rms_relative']:.4f}  "
              f"a_y0={record['fit_params']['a_y0']:.4g} "
              f"y_sep={record['fit_params']['y_sep']:.4g} "
              f"a_y1={record['fit_params']['a_y1']:.4g}")

    base_top = {v: ped_top(x0, np.asarray(getattr(phys0, v).values))
                for v in variables}
    print(f"[test] nominal pedestal-top values: "
          + ", ".join(f"{v}={base_top[v]:.4g}" for v in variables))

    rows = []
    for var in variables:
        for s in scales:
            tag = f"{var}_{s:.3f}"
            savedir = os.path.join(workdir, tag)
            print(f"\n=== {tag} " + "=" * 50)

            phys = phys0.apply_mtanh_ped(var, fit=fits[var], scale_height=s,
                                         enforce_quasineutrality=True, qz=6.0)

            # Q1 is answered before cheaseBS even runs: the mapping is a property
            # of the transform, not of the equilibrium solve.
            x = phys.rhot.values if hasattr(phys.rhot, "values") else np.asarray(phys.rhot)
            measured = {v: ped_top(x, np.asarray(getattr(phys, v).values))
                        for v in ("ne", "Te", "ni", "nz") if hasattr(phys, v)}
            top_ratio = measured[var] / base_top[var]
            print(f"[test] scale {s:.3f} -> pedestal-top ratio {top_ratio:.4f} "
                  f"({(top_ratio - 1) * 100:+.2f}%)")

            row = {"var": var, "scale": s, "tag": tag,
                   "ped_top": measured, "ped_top_ratio": top_ratio,
                   "savedir": savedir}

            try:
                phys.output_gfile(
                    savedir=savedir, run_cheasebs=True,
                    cheasebs_acceptance=CheasebsAcceptance.production(
                        analysis_radii=GENE_RADII),
                    cheasebs_strict=False, comment=tag, **CHEASEBS_OVERRIDES)
                rec = read_acceptance(savedir) or {}
            except Exception as exc:
                row["error"] = f"{type(exc).__name__}: {exc}"
                print(f"[test] FAILED: {row['error']}")
                rows.append(row)
                continue

            summary_path = os.path.join(savedir, "convergence_summary.json")
            summary = {}
            if os.path.exists(summary_path):
                summary = json.load(open(summary_path))

            row.update({
                "accepted": rec.get("accepted"),
                "ip_error_rel": rec.get("ip_error_rel"),
                "q_errors_rel": rec.get("q_errors_rel"),
                "iterations": summary.get("iterations"),
                "converged": summary.get("final_converged"),
                "bootstrap_change_rel": summary.get("final_bootstrap_change_rel"),
                "q_change_rel": summary.get("final_q_change_rel"),
                "final_ip_a": summary.get("final_ip_a"),
                "target_ip_a": summary.get("target_ip_a"),
            })
            rows.append(row)

    # ---- report -------------------------------------------------------
    out = os.path.join(workdir, "scale_symmetry_results.json")
    with open(out, "w") as f:
        json.dump({"dirpath": args.dirpath, "scales": scales,
                   "variables": variables, "nominal_fit": nominal,
                   "nominal_ped_top": base_top,
                   "cheasebs_overrides": CHEASEBS_OVERRIDES,
                   "rows": rows}, f, indent=1, default=str)

    print("\n" + "=" * 78)
    print("Q1 — scale factor to pedestal-top mapping")
    head = f"{'var':<4} {'scale':>7} {'top ratio':>10} {'top %':>8} {'slope':>8}"
    print(head); print("-" * len(head))
    for r in rows:
        if "ped_top_ratio" not in r:
            continue
        d = r["scale"] - 1.0
        slope = (r["ped_top_ratio"] - 1.0) / d if abs(d) > 1e-12 else float("nan")
        print(f"{r['var']:<4} {r['scale']:>7.3f} {r['ped_top_ratio']:>10.4f} "
              f"{(r['ped_top_ratio']-1)*100:>7.2f}% {slope:>8.4f}")
    print("  A constant slope column means the mapping is linear and one number "
          "per variable converts the survey's absolute bounds into scale factors.")

    print("\nQ2 — response symmetry")
    head = (f"{'var':<4} {'scale':>7} {'iters':>6} {'conv':>6} {'Ip err':>9} "
            f"{'final Ip':>12} {'bs_change':>11}")
    print(head); print("-" * len(head))
    for r in rows:
        if "iterations" not in r:
            print(f"{r['var']:<4} {r['scale']:>7.3f}   {r.get('error','no result')}")
            continue
        ip = r["ip_error_rel"]
        print(f"{r['var']:<4} {r['scale']:>7.3f} {str(r['iterations']):>6} "
              f"{str(r['converged']):>6} "
              f"{(f'{ip:.4%}' if isinstance(ip,(int,float)) else 'n/a'):>9} "
              f"{(r['final_ip_a'] or float('nan')):>12.1f} "
              f"{str(r['bootstrap_change_rel']):>11}")

    print("\n  Compare each +delta row against its -delta partner. Equal iteration")
    print("  counts and equal-and-opposite Ip shifts mean the response is")
    print("  symmetric and the 19-vs-2 split was specific to that one point.")
    print("  A systematic difference points at a one-sided mechanism -- the")
    print("  floored p_fast split, nz clipping under quasineutrality at low")
    print("  density, or qspec pinning the current.")
    print(f"\n[test] results: {out}")
    print(f"[test] per-point iteration_errors.png under {workdir}/<tag>/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
