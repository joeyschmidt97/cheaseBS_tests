"""Package test U (132543 NSTX-U projection, B0 0.439 -> 0.9 T) for outside review.

Two jobs:

1. Preliminary check. Reads the source EFIT, the final cheaseBS EQDSK, the
   scaled and reference kinetic profiles, and (when present) the unscaled
   test-K baseline, then checks the scaling the run was meant to impose:

       B0, Ip, F_edge          x lambda           (lambda = 0.9 / 0.439262)
       Ip / B0                 held               (so q should be held)
       T_e, T_i, T_z           T_sep + lambda^2 (T - T_sep)
       n_e, n_i, n_z           held
       p_total                 p_thermal(new) + p_fast(source, floored at 0)
       beta_T, beta_N, beta_p  x s_p / lambda^2   (s_p = <p>_V ratio)
       q profile, l_i, R_axis, boundary shape     ideally held

   Every check is PASS / WARN / FAIL with its tolerance written next to it.
   Results go to review/ beside this script (check_report.md, checks.json,
   PNGs) so they can be committed with the run.

2. Tarball. Stages the final case files, the iteration_00 and last-iteration
   directories from the NERSC scratch run directory, an error digest over
   every CHEASE log (plus the logs of crashed attempts named in console.log),
   and the review/ outputs, then writes <out>/<name>.tar.gz.

Run on NERSC from the cheaseBS_tests repo root (the scratch run directory is
only there; locally the iteration copy is skipped with a warning):

    python data/cheasebs_132xxx/U_nstxu_0.9T/pack_for_review.py
    python data/cheasebs_132xxx/U_nstxu_0.9T/pack_for_review.py --analysis-only
    python data/cheasebs_132xxx/U_nstxu_0.9T/pack_for_review.py --skip-ogyropsi --out $SCRATCH

Needs numpy and matplotlib only (TPED conda env is fine).
"""

import argparse
import csv
import datetime
import json
import re
import shutil
import sys
import tarfile
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath

HERE = Path(__file__).resolve().parent
MU0 = 4e-7 * np.pi
QE = 1.602176634e-19

DEFAULT_CASE = HERE / "132543" / "T_scale_4.198_sep"
DEFAULT_BASELINE = (HERE.parent / "K_rhop_tensbnd" / "132543" / "Te_ped_scale_1.000"
                    / "EQDSK_COCOS_02_POS_SOURCE_SIGNS.OUT")

# CHEASE log lines worth a reviewer's attention. "chease_input.xml file
# missing" is printed on every run and is benign, so it is counted apart.
FATAL_PATTERNS = [r"BRACKETING NOT SUCCEEDED", r"error after", r"cannot find",
                  r"Segmentation", r"SIGSEGV", r"forrtl", r"\bNaN\b", r"STOP\b",
                  r"did not produce", r"Traceback", r"Error\b", r"crashed"]
WARN_PATTERNS = [r"WARNING"]
BENIGN_PATTERNS = [r"chease_input\.xml"]


# --------------------------------------------------------------------------
# readers
# --------------------------------------------------------------------------

_NUM = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?")


def read_geqdsk(path):
    """Minimal G-EQDSK reader. Handles fields that run together ('1.0E+00-2.0E-01')."""
    text = Path(path).read_text()
    head, rest = text.split("\n", 1)
    nw, nh = (int(t) for t in head.split()[-2:])
    v = [float(t) for t in _NUM.findall(rest)]
    g = dict(nw=nw, nh=nh, header=head.strip(), path=str(path))
    (g["rdim"], g["zdim"], g["rcentr"], g["rleft"], g["zmid"],
     g["rmaxis"], g["zmaxis"], g["simag"], g["sibry"], g["bcentr"],
     g["current"]) = v[:11]
    i = 20
    for key in ("fpol", "pres", "ffprim", "pprime"):
        g[key] = np.array(v[i:i + nw]); i += nw
    g["psirz"] = np.array(v[i:i + nw * nh]).reshape(nh, nw); i += nw * nh
    g["qpsi"] = np.array(v[i:i + nw]); i += nw
    nbbbs, limitr = int(v[i]), int(v[i + 1]); i += 2
    bb = np.array(v[i:i + 2 * nbbbs]).reshape(nbbbs, 2)
    g["rbbbs"], g["zbbbs"] = bb[:, 0], bb[:, 1]
    g["psin_1d"] = np.linspace(0.0, 1.0, nw)
    g["R"] = g["rleft"] + np.linspace(0.0, g["rdim"], nw)
    g["Z"] = g["zmid"] - g["zdim"] / 2 + np.linspace(0.0, g["zdim"], nh)
    return g


def read_profile(path):
    """GENE profiles_* file: rho_tor, rho_pol, T [keV], n [1e19 m^-3]."""
    a = np.loadtxt(path, comments="#")
    return dict(rhot=a[:, 0], rhop=a[:, 1], T=a[:, 2], n=a[:, 3])


def load_json(path):
    try:
        return json.loads(Path(path).read_text())
    except (OSError, ValueError):
        return {}


# --------------------------------------------------------------------------
# equilibrium scalars
# --------------------------------------------------------------------------

def scalars(g):
    """Global quantities from one EQDSK. Definitions match the 2026-10-24 spot check."""
    R, Z = np.meshgrid(g["R"], g["Z"])
    dR, dZ = g["R"][1] - g["R"][0], g["Z"][1] - g["Z"][0]
    psin = (g["psirz"] - g["simag"]) / (g["sibry"] - g["simag"])
    bnd = np.column_stack([g["rbbbs"], g["zbbbs"]])
    inside = MplPath(bnd).contains_points(np.column_stack([R.ravel(), Z.ravel()]))
    inside = inside.reshape(R.shape) & (psin >= 0) & (psin <= 1)

    x = g["psin_1d"]
    pp = np.interp(psin, x, g["pprime"])
    ffp = np.interp(psin, x, g["ffprim"])
    p2d = np.interp(psin, x, g["pres"])
    jphi = R * pp + ffp / (MU0 * R)
    ip_int = abs(np.sum(jphi[inside]) * dR * dZ)

    dpsi_dZ, dpsi_dR = np.gradient(g["psirz"], dZ, dR)
    bp2 = (dpsi_dR ** 2 + dpsi_dZ ** 2) / R ** 2
    dV = 2 * np.pi * R * dR * dZ
    vol = np.sum(dV[inside])
    p_vol = np.sum((p2d * dV)[inside]) / vol

    rb, zb = g["rbbbs"], g["zbbbs"]
    rmax, rmin, zmax, zmin = rb.max(), rb.min(), zb.max(), zb.min()
    a = (rmax - rmin) / 2
    rgeo = (rmax + rmin) / 2
    kappa = (zmax - zmin) / (2 * a)
    d_up = (rgeo - rb[np.argmax(zb)]) / a
    d_lo = (rgeo - rb[np.argmin(zb)]) / a
    rc, zc = np.append(rb, rb[0]), np.append(zb, zb[0])
    L = np.sum(np.hypot(np.diff(rc), np.diff(zc)))

    ip = abs(g["current"])
    b0 = abs(g["bcentr"])
    bpa = MU0 * ip / L
    beta_t = 2 * MU0 * p_vol / b0 ** 2
    q = np.abs(g["qpsi"])
    imin = int(np.argmin(q))
    return dict(
        B0_T=b0, Ip_MA=ip / 1e6, Ip_int_MA=ip_int / 1e6,
        F_edge_Tm=abs(g["fpol"][-1]), psi_span_Wb=abs(g["sibry"] - g["simag"]),
        p_vol_kPa=p_vol / 1e3, p0_kPa=g["pres"][0] / 1e3,
        beta_T_pct=100 * beta_t,
        beta_N=100 * beta_t * a * b0 / (ip / 1e6),
        beta_p=2 * MU0 * p_vol / bpa ** 2,
        li3=2 * np.sum((bp2 * dV)[inside]) / (MU0 ** 2 * ip ** 2 * rgeo),
        q0=q[0], qmin=q[imin], psin_qmin=x[imin],
        q95=float(np.interp(0.95, x, q)),
        R_axis_m=g["rmaxis"], Z_axis_m=g["zmaxis"],
        a_m=a, R_geo_m=rgeo, kappa=kappa, delta_up=d_up, delta_lo=d_lo,
        volume_m3=vol, L_lcfs_m=L,
    )


# --------------------------------------------------------------------------
# checks
# --------------------------------------------------------------------------

class Checks:
    def __init__(self):
        self.rows = []

    def add(self, name, value, expected, tol_pass, tol_warn=None, note="", rel=True):
        """Status from |value/expected - 1| (rel) or |value - expected| (abs)."""
        if value is None or expected is None:
            dev, status = None, "SKIP"
        else:
            dev = abs(value / expected - 1) if rel else abs(value - expected)
            if dev <= tol_pass:
                status = "PASS"
            elif tol_warn is not None and dev <= tol_warn:
                status = "WARN"
            else:
                status = "FAIL" if tol_warn is None else "WARN"
        self.rows.append(dict(check=name, value=value, expected=expected, deviation=dev,
                              rel=rel, tol_pass=tol_pass, tol_warn=tol_warn,
                              status=status, note=note))

    def flag(self, name, ok, note="", warn_only=False):
        status = "PASS" if ok else ("WARN" if warn_only else "FAIL")
        self.rows.append(dict(check=name, value=None, expected=None, deviation=None,
                              rel=True, tol_pass=None, tol_warn=None,
                              status=status, note=note))

    def counts(self):
        out = {}
        for r in self.rows:
            out[r["status"]] = out.get(r["status"], 0) + 1
        return out


def _fmt(x):
    if x is None:
        return "-"
    if isinstance(x, float):
        return f"{x:.4g}"
    return str(x)


def run_analysis(case, source_g, baseline_g, run_dir, outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    tlog = load_json(case / "TPED_transform_log.json")
    cfg = load_json(case / "cheasebs_run_config.json")
    acc = load_json(case / "cheasebs_acceptance.json")
    conv = load_json(case / "convergence_summary.json")

    final_g_path = next(iter(sorted(case.glob("*.geqdsk"))), None) \
        or next(iter(sorted(case.glob("EQDSK*.OUT"))), None)
    if final_g_path is None:
        sys.exit(f"no final EQDSK (*.geqdsk / EQDSK*.OUT) in {case}")

    src = read_geqdsk(source_g)
    fin = read_geqdsk(final_g_path)
    base = read_geqdsk(baseline_g) if baseline_g and Path(baseline_g).exists() else None
    S = {"source": scalars(src), "final": scalars(fin)}
    if base is not None:
        S["baseline_K"] = scalars(base)

    # Scale factors from the run itself, not hard-coded.
    lam = abs(cfg.get("target_bt_t", -0.9)) / S["source"]["B0_T"]
    lam_ip = cfg.get("target_ip_a", lam * S["source"]["Ip_MA"] * 1e6) / (S["source"]["Ip_MA"] * 1e6)
    th = next((t for t in tlog.get("transform_history", []) if t.get("transform") == "apply_t_scale"), {})
    t_factor = th.get("factor", lam ** 2)
    t_sep_ev = th.get("t_sep", {})

    C = Checks()
    sF, sS = S["final"], S["source"]
    # Ratio checks run against the unscaled cheaseBS baseline when it exists: it
    # shares CHEASE's boundary (EFIT's lower X-point is rounded), so the ratio
    # isolates the scaling. Falls back to the source EFIT.
    rl = "baseline K" if "baseline_K" in S else "source"
    sR = S.get("baseline_K", sS)
    note_lam = f"lambda = {lam:.4f}"
    C.add("B0 ratio final/source", sF["B0_T"] / sS["B0_T"], lam, 2e-3, note=note_lam)
    C.add("Ip ratio final/source", sF["Ip_MA"] / sS["Ip_MA"], lam_ip, 2e-3, 5e-3, note="target tol_ip_rel = 1e-3")
    C.add("Ip/B0 held (final/source)", (sF["Ip_MA"] / sF["B0_T"]) / (sS["Ip_MA"] / sS["B0_T"]), 1.0, 2e-3, 5e-3)
    C.add("Ip vs cheaseBS target", sF["Ip_MA"] * 1e6, cfg.get("target_ip_a"), 1e-3, 5e-3)
    C.add("integral j_phi dA vs header Ip", sF["Ip_int_MA"], sF["Ip_MA"], 1e-2, 3e-2,
          note="grid integral of R p' + FF'/(mu0 R) inside LCFS")
    C.add(f"F_edge ratio final/{rl}", sF["F_edge_Tm"] / sR["F_edge_Tm"], lam, 1e-2, 2e-2)
    C.add(f"psi span ratio final/{rl}", sF["psi_span_Wb"] / sR["psi_span_Wb"], lam, 3e-2, 6e-2,
          note="exact similarity gives lambda; current-profile change moves it")

    # ---- kinetic profiles: T_new = T_sep + f (T_ref - T_sep), n held
    prof = {}
    ref_dir = case / "reference_profiles"
    for sp, key in (("e", "Te"), ("i", "Ti"), ("z", "Tz")):
        new_p, ref_p = case / f"profiles_{sp}", ref_dir / f"REF_profiles_{sp}"
        if not (new_p.exists() and ref_p.exists()):
            C.flag(f"profiles_{sp} present", False, f"missing {new_p.name} or {ref_p}")
            continue
        new, ref = read_profile(new_p), read_profile(ref_p)
        tsep = t_sep_ev.get(key, 1e3 * ref["T"][-1]) / 1e3
        t_ref = np.interp(new["rhop"], ref["rhop"], ref["T"])
        n_ref = np.interp(new["rhop"], ref["rhop"], ref["n"])
        t_exp = tsep + t_factor * (t_ref - tsep)
        prof[sp] = dict(new=new, ref=ref, t_exp=t_exp, tsep=tsep)
        C.add(f"{key} = T_sep + {t_factor:.3f}(T - T_sep), max rel dev",
              float(np.max(np.abs(new["T"] / t_exp - 1))), 0.0, 1e-3, 1e-2, rel=False,
              note=f"T_sep = {1e3 * tsep:.1f} eV")
        C.add(f"n_{sp} held, max rel dev", float(np.max(np.abs(new["n"] / n_ref - 1))), 0.0,
              1e-6, 1e-3, rel=False)
        C.add(f"{key} separatrix held", new["T"][-1], ref["T"][-1], 1e-3, 1e-2)

    # ---- pressure bookkeeping: p_final = p_th(new) + p_fast(source)
    x = src["psin_1d"]
    if len(prof) == 3:
        def p_th(key):  # Pa on the source psi_N grid
            tot = 0.0
            for sp in ("e", "i", "z"):
                d = prof[sp]
                t = d["new"]["T"] if key == "new" else d["ref"]["T"]
                rp = d["new"]["rhop"] if key == "new" else d["ref"]["rhop"]
                n = d["new"]["n"] if key == "new" else d["ref"]["n"]
                tot = tot + np.interp(np.sqrt(x), rp, n * t) * 1e19 * 1e3 * QE
            return tot
        p_th_ref, p_th_new = p_th("ref"), p_th("new")
        p_fast = np.clip(src["pres"] - p_th_ref, 0.0, None)
        p_exp = p_th_new + p_fast
        p_fin = np.interp(x, fin["psin_1d"], fin["pres"])
        mask = x <= 0.98
        dev = float(np.max(np.abs(p_fin[mask] - p_exp[mask])) / p_exp[0])
        C.add("p_final = p_th(new) + p_fast(src), max |dp|/p_exp(0), psiN<=0.98",
              dev, 0.0, 2e-2, 5e-2, rel=False,
              note=f"p_fast/p on axis = {p_fast[0] / src['pres'][0]:.2f}")
        prof["_p"] = dict(x=x, p_th_ref=p_th_ref, p_th_new=p_th_new, p_fast=p_fast,
                          p_exp=p_exp, p_fin=p_fin)

    s_p = sF["p_vol_kPa"] / sR["p_vol_kPa"]
    exp_beta = s_p / lam ** 2
    C.add(f"beta_T ratio final/{rl} vs s_p/lambda^2", sF["beta_T_pct"] / sR["beta_T_pct"], exp_beta, 1e-2, 3e-2,
          note=f"s_p = <p>_V ratio = {s_p:.3f}; exact-beta match would be 1")
    C.add(f"beta_N ratio final/{rl} vs s_p/(s_B s_I)", sF["beta_N"] / sR["beta_N"], s_p / (lam * lam_ip), 1e-2, 3e-2)
    C.add(f"beta_p ratio final/{rl} vs s_p/s_I^2", sF["beta_p"] / sR["beta_p"], s_p / lam_ip ** 2, 1e-2, 3e-2)
    C.add("beta_T held (thermal-only scaling target)", sF["beta_T_pct"] / sR["beta_T_pct"], 1.0, 2e-2, 0.25,
          note="WARN expected: p_fast held at source value")

    # ---- q and shape (ideally held under Ip/B0 fixed)
    C.add(f"q95 final vs {rl}", sF["q95"], sR["q95"], 2e-2, 5e-2)
    C.add(f"q0 final vs {rl}", sF["q0"], sR["q0"], 2e-2, 0.15,
          note="bootstrap recomputed at lower nu*, so q0 may move")
    C.add(f"qmin final vs {rl}", sF["qmin"], sR["qmin"], 2e-2, 0.15, note=f"qmin at psiN = {sF['psin_qmin']:.3f}")
    ref_g = base if base is not None else src
    qf = np.interp(x, fin["psin_1d"], np.abs(fin["qpsi"]))
    qr = np.interp(x, ref_g["psin_1d"], np.abs(ref_g["qpsi"]))
    band = (x >= 0.1) & (x <= 0.9)
    C.add(f"max |q_final/q_{rl} - 1|, psiN 0.1-0.9", float(np.max(np.abs(qf[band] / qr[band] - 1))),
          0.0, 2e-2, 0.10, rel=False)
    C.add(f"l_i(3) final vs {rl}", sF["li3"], sR["li3"], 3e-2, 0.10)
    C.add(f"R_axis final vs {rl} [m]", sF["R_axis_m"], sR["R_axis_m"], 5e-3, 2e-2, rel=False)
    C.add(f"minor radius a vs {rl}", sF["a_m"], sR["a_m"], 5e-3, 2e-2)
    C.add(f"elongation kappa vs {rl}", sF["kappa"], sR["kappa"], 1e-2, 3e-2)
    C.add("q95 final vs source EFIT", sF["q95"], sS["q95"], 2e-2, 5e-2)

    # ---- convergence and cheaseBS's own acceptance
    C.flag("cheaseBS converged", bool(conv.get("final_converged")),
           f"{conv.get('iterations')} iterations; final Ip err {conv.get('final_ip_error_rel', float('nan')):.2e}, "
           f"bs change {conv.get('final_bootstrap_change_rel', float('nan')):.2e}, "
           f"q change {conv.get('final_q_change_rel', float('nan')):.2e}")
    C.flag("cheaseBS acceptance", bool(acc.get("accepted")),
           f"q95 err {acc.get('q95_error_rel', float('nan')):.3f}, q-edge err "
           f"{acc.get('q_edge_error_rel', float('nan')):.3f} (vs cheaseBS internal q target)")

    # ---- iteration 00 vs last, from the scratch run directory when present
    iters = {}
    hist = []
    log_csv = case / "iteration_log.csv"
    if log_csv.exists():
        with open(log_csv) as fh:
            hist = list(csv.DictReader(fh))
    if run_dir is not None and run_dir.exists():
        last = last_iteration(run_dir, conv)
        for name in ("iteration_00", last):
            p = next(iter(sorted((run_dir / name / "artifacts").glob("EQDSK*SOURCE_SIGNS*.OUT"))), None) \
                or next(iter(sorted((run_dir / name / "artifacts").glob("EQDSK*.OUT"))), None)
            if p is not None:
                iters[name] = scalars(read_geqdsk(p))
        if "iteration_00" in iters:
            S["iteration_00"] = iters["iteration_00"]

    write_plots(outdir, src, fin, base, prof, hist, lam, t_factor, S)
    write_report(outdir, C, S, lam, lam_ip, t_factor, s_p, case, final_g_path, source_g, baseline_g, run_dir, iters)
    return C


def last_iteration(run_dir, conv):
    if "final_iteration" in conv:
        name = f"iteration_{int(conv['final_iteration']):02d}"
        if (run_dir / name).exists():
            return name
    dirs = sorted(p.name for p in run_dir.glob("iteration_[0-9]*") if p.is_dir())
    return dirs[-1] if dirs else None


# --------------------------------------------------------------------------
# outputs
# --------------------------------------------------------------------------

def write_plots(outdir, src, fin, base, prof, hist, lam, t_factor, S):
    x = src["psin_1d"]
    xf = fin["psin_1d"]
    fig, ax = plt.subplots(2, 3, figsize=(15, 8.5))

    a = ax[0, 0]
    a.plot(x, np.abs(src["qpsi"]), "k", label="source EFIT")
    if base is not None:
        a.plot(base["psin_1d"], np.abs(base["qpsi"]), "C0--", label="baseline K (unscaled)")
    a.plot(xf, np.abs(fin["qpsi"]), "C3", label="final U")
    a.set(xlabel=r"$\psi_N$", ylabel="q", title="q profile (Ip/B0 held: ideally unchanged)", ylim=(0, 15))
    a.legend()

    a = ax[0, 1]
    if "_p" in prof:
        P = prof["_p"]
        a.plot(x, src["pres"] / 1e3, "k", label="source EFIT p")
        a.plot(x, P["p_fast"] / 1e3, "k:", label="p_fast (held)")
        a.plot(x, P["p_exp"] / 1e3, "C2--", lw=2, label="expected p_th(new)+p_fast")
        a.plot(x, P["p_fin"] / 1e3, "C3", label="final U")
    a.set(xlabel=r"$\psi_N$", ylabel="p [kPa]", title="pressure bookkeeping")
    a.legend()

    # Against the CHEASE baseline when present: EFIT's FF' crosses zero, so a
    # ratio to it blows up without saying anything about the scaling.
    ref, rl = (base, "K") if base is not None else (src, "src")
    a = ax[0, 2]
    pf = np.interp(x, xf, fin["pprime"])
    ff = np.interp(x, xf, fin["ffprim"])
    pr = np.interp(x, ref["psin_1d"], ref["pprime"])
    fr = np.interp(x, ref["psin_1d"], ref["ffprim"])
    with np.errstate(divide="ignore", invalid="ignore"):
        a.plot(x, pf / pr / lam, "C1", label=rf"$p'_U / p'_{{{rl}}} / \lambda$")
        a.plot(x, ff / fr / lam, "C4", label=rf"$FF'_U / FF'_{{{rl}}} / \lambda$")
    a.axhline(1, color="k", lw=0.5)
    a.set(xlabel=r"$\psi_N$", ylabel="ratio", ylim=(-1, 4),
          title="source functions / exact-similarity scaling (=1)")
    a.legend()

    a = ax[1, 0]
    for sp, c in (("e", "C0"), ("i", "C3"), ("z", "C2")):
        if sp in prof:
            d = prof[sp]
            a.plot(d["ref"]["rhop"], d["ref"]["T"], c + ":", label=f"T{sp} ref")
            a.plot(d["new"]["rhop"], d["new"]["T"], c, label=f"T{sp} U")
            a.plot(d["new"]["rhop"], d["t_exp"], "k--", lw=0.7)
    a.set(xlabel=r"$\rho_{pol}$", ylabel="T [keV]",
          title=f"T = T_sep + {t_factor:.3f}(T - T_sep)  (dashed: expected)")
    a.legend(fontsize=8)

    a = ax[1, 1]
    a.plot(src["rbbbs"], src["zbbbs"], "k", label="source EFIT LCFS")
    a.plot(fin["rbbbs"], fin["zbbbs"], "C3--", label="final U LCFS")
    a.plot(src["rmaxis"], src["zmaxis"], "k+", ms=10)
    a.plot(fin["rmaxis"], fin["zmaxis"], "C3x", ms=10)
    a.set(xlabel="R [m]", ylabel="Z [m]", title="boundary and magnetic axis", aspect="equal")
    a.legend(fontsize=8)

    a = ax[1, 2]
    if hist:
        it = [int(h["iteration"]) for h in hist]
        for key, c in (("ip_error_rel", "C0"), ("bootstrap_change_rel", "C1"), ("q_change_rel", "C2")):
            vals = [float(h[key]) if h[key] else np.nan for h in hist]
            a.semilogy(it, vals, "o-", color=c, label=key)
        a.axhline(1e-3, color="k", ls=":", label="tolerance 1e-3")
        a2 = a.twinx()
        a2.plot(it, [float(h["amplitude_in"]) for h in hist], "k^-", alpha=0.4)
        a2.set_ylabel("driven-current amplitude")
        a.legend(fontsize=8)
    a.set(xlabel="iteration", title="cheaseBS convergence history")

    fig.suptitle(f"Test U 132543: NSTX-U projection, lambda = {lam:.3f}")
    fig.tight_layout()
    fig.savefig(outdir / "scaling_check.png", dpi=130)
    plt.close(fig)


def write_report(outdir, C, S, lam, lam_ip, t_factor, s_p, case, final_g, source_g, baseline_g, run_dir, iters):
    counts = C.counts()
    cols = list(S)
    keys = list(S["source"])
    L = [f"# Test U scaling check: 132543 at 0.9 T",
         "",
         f"Generated {datetime.datetime.now():%Y-%m-%d %H:%M} by `pack_for_review.py`.",
         "",
         f"- Final EQDSK: `{final_g}`",
         f"- Source EFIT: `{source_g}`",
         f"- Baseline K (unscaled, same settings): `{baseline_g if baseline_g and Path(baseline_g).exists() else 'not found'}`",
         f"- Scratch run dir: `{run_dir if run_dir and run_dir.exists() else 'not available here'}`",
         f"- lambda (B0) = {lam:.4f}, lambda (Ip) = {lam_ip:.4f}, T factor = {t_factor:.4f}, "
         f"s_p = <p>_V ratio = {s_p:.3f}",
         "",
         "**Summary:** " + ", ".join(f"{v} {k}" for k, v in sorted(counts.items())),
         "",
         "## Checks",
         "",
         "| status | check | value | expected | deviation | pass / warn tol | note |",
         "|---|---|---|---|---|---|---|"]
    for r in C.rows:
        tol = "-" if r["tol_pass"] is None else f"{r['tol_pass']:g} / {_fmt(r['tol_warn'])}"
        L.append(f"| {r['status']} | {r['check']} | {_fmt(r['value'])} | {_fmt(r['expected'])} | "
                 f"{_fmt(r['deviation'])} | {tol} | {r['note']} |")
    L += ["", "## Equilibrium scalars", "",
          "| quantity | " + " | ".join(cols) + " | final/source |",
          "|---|" + "---|" * (len(cols) + 1)]
    for k in keys:
        ratio = S["final"][k] / S["source"][k] if S["source"][k] else float("nan")
        L.append(f"| {k} | " + " | ".join(_fmt(S[c][k]) for c in cols) + f" | {ratio:.4f} |")
    if iters:
        L += ["", "## Iteration 00 vs last iteration", "",
              "| quantity | " + " | ".join(iters) + " |", "|---|" + "---|" * len(iters)]
        for k in ("Ip_MA", "q0", "qmin", "q95", "li3", "beta_N", "R_axis_m"):
            L.append(f"| {k} | " + " | ".join(_fmt(v[k]) for v in iters.values()) + " |")
    L += ["", "## Definitions", "",
          "- beta_T = 2 mu0 <p>_V / B0^2 with total EQDSK pressure (thermal + fast).",
          "- beta_N = beta_T[%] a B0 / Ip[MA], a = (Rmax - Rmin)/2 of the LCFS.",
          "- beta_p uses B_p,a = mu0 Ip / L_LCFS.",
          "- l_i(3) = 2 int B_p^2 dV / (mu0^2 Ip^2 R_geo).",
          "- p_fast = max(p_EFIT - p_thermal(reference profiles), 0), held at source value by cheaseBS.",
          "- WARN on beta_T held and on q0/l_i is expected for this case; see the run notes.",
          "", "![scaling check](scaling_check.png)", ""]
    (outdir / "check_report.md").write_text("\n".join(L), encoding="utf-8")
    (outdir / "checks.json").write_text(json.dumps(
        dict(lambda_B0=lam, lambda_Ip=lam_ip, t_factor=t_factor, s_p=s_p,
             counts=counts, checks=C.rows, scalars=S, iterations=iters), indent=2, default=float))


# --------------------------------------------------------------------------
# error digest
# --------------------------------------------------------------------------

def scan_log(path):
    hits = {"fatal": [], "warning": [], "benign": 0}
    try:
        lines = Path(path).read_text(errors="replace").splitlines()
    except OSError:
        return hits
    for i, line in enumerate(lines, 1):
        if any(re.search(p, line) for p in BENIGN_PATTERNS):
            hits["benign"] += 1
        elif any(re.search(p, line) for p in FATAL_PATTERNS):
            hits["fatal"].append((i, line.strip()))
        elif any(re.search(p, line) for p in WARN_PATTERNS):
            hits["warning"].append((i, line.strip()))
    return hits


def error_digest(case, test_dir, run_dir, errdir):
    """Scan every CHEASE log of the run, plus logs of crashed attempts listed in console.log."""
    errdir.mkdir(parents=True, exist_ok=True)
    logs = []
    if run_dir is not None and run_dir.exists():
        logs += sorted(run_dir.rglob("*.log"))
    logs += sorted((case / "chease_logs").glob("*.log"))

    console = test_dir / "console.log"
    other_runs = []
    if console.exists():
        text = console.read_text(errors="replace")
        for m in sorted(set(re.findall(r"(/pscratch/\S+/cheaseBS_runs/\d+/[^\s/()]+)", text))):
            p = Path(m)
            if run_dir is not None and p.resolve() == run_dir.resolve():
                continue
            other_runs.append(p)
            if p.exists():
                dest = errdir / "other_attempts" / p.name
                for lg in p.rglob("*.log"):
                    tgt = dest / lg.relative_to(p)
                    tgt.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(lg, tgt)
                    logs.append(lg)
        shutil.copy2(console, errdir / "console.log")

    L = ["# Error digest", "",
         f"{len(logs)} log files scanned. 'chease_input.xml file missing' is printed by every "
         "CHEASE run and is counted as benign.", ""]
    tb = []
    if console.exists():
        lines = console.read_text(errors="replace").splitlines()
        for i, line in enumerate(lines):
            if "Traceback" in line or "crashed" in line or "CHEASE did not produce" in line:
                tb.append(f"console.log:{i + 1}: {line.strip()}")
    if tb:
        L += ["## console.log failures (earlier attempts of this test)", ""] + [f"- `{t}`" for t in tb] + [""]
    if other_runs:
        L += ["## Other run directories named in console.log", ""]
        L += [f"- `{p}` ({'logs copied to other_attempts/' if p.exists() else 'not found here'})"
              for p in other_runs] + [""]
    L += ["## Per-log findings", "", "| log | fatal | warning | benign |", "|---|---|---|---|"]
    detail = []
    for lg in logs:
        h = scan_log(lg)
        L.append(f"| `{lg}` | {len(h['fatal'])} | {len(h['warning'])} | {h['benign']} |")
        if h["fatal"] or h["warning"]:
            detail.append(f"### `{lg}`")
            detail += [f"- FATAL l.{i}: `{s}`" for i, s in h["fatal"][:20]]
            detail += [f"- warn l.{i}: `{s}`" for i, s in h["warning"][:20]]
            detail.append("")
    L += ["", "## Matched lines", ""] + detail
    (errdir / "error_digest.md").write_text("\n".join(L), encoding="utf-8")
    return tb


# --------------------------------------------------------------------------
# staging and tar
# --------------------------------------------------------------------------

def stage(case, test_dir, run_dir, review_dir, stage_dir, skip_ogyropsi, with_baseline, conv):
    ignore = shutil.ignore_patterns("ogyropsi.dat") if skip_ogyropsi else None
    shutil.copytree(case, stage_dir / "final" / case.name, ignore=ignore)
    tl = stage_dir / "test_level"
    tl.mkdir(parents=True)
    for f in test_dir.iterdir():
        if f.is_file() and f.suffix != ".py":
            shutil.copy2(f, tl / f.name)

    copied = []
    if run_dir is not None and run_dir.exists():
        last = last_iteration(run_dir, conv)
        for name in dict.fromkeys(["iteration_00", last]):
            if name and (run_dir / name).exists():
                shutil.copytree(run_dir / name, stage_dir / "iterations" / name, ignore=ignore)
                copied.append(name)
        top = stage_dir / "run_dir_top"
        top.mkdir()
        for f in run_dir.iterdir():
            if f.is_file() and f.stat().st_size < 5e6:
                shutil.copy2(f, top / f.name)
        if with_baseline and (run_dir / "cheasebs_baseline").exists():
            shutil.copytree(run_dir / "cheasebs_baseline", stage_dir / "cheasebs_baseline", ignore=ignore)
    else:
        print(f"WARNING: scratch run dir not found ({run_dir}); iteration directories not packed")

    shutil.copytree(review_dir, stage_dir / "review")
    return copied


def write_readme(stage_dir, copied, run_dir, counts):
    txt = f"""# Test U review package: NSTX 132543 projected to NSTX-U (0.9 T)

Fixed-boundary CHEASE (cheaseBS iterative bootstrap closure) of 132543 at 700 ms,
scaled to B0 = 0.9 T. Ip and B0 both x lambda = 2.049 (Ip/B0 held, Ip = 1.94 MA).
Te, Ti, Tz -> T_sep + lambda^2 (T - T_sep) with separatrix values and density held.
Fast-ion pressure held at the source value. Settings: jparallel replay on rho_pol,
NCSCAL = 4, TENSBND = -1.0, bootstrap_mix 0.35, initial amplitude = lambda.

Start here
- review/check_report.md   preliminary scaling checks ({', '.join(f'{v} {k}' for k, v in sorted(counts.items()))})
- review/scaling_check.png q, pressure, source functions, T profiles, boundary, convergence
- errors/error_digest.md   every CHEASE log scanned, plus the crashed first attempts

Layout
- final/<case>/            final EQDSK (*.geqdsk), profiles_e/i/z (scaled), reference_profiles/
                           (unscaled), profiles.iterdb, run config, acceptance, convergence,
                           iteration_log.csv, TPED transform log, chease_logs/
- test_level/              namelist, test summary.json, console.log, comparison figures
- iterations/              {', '.join(copied) if copied else 'not included (run dir unavailable)'}
                           full CHEASE artifacts (EQDSK, EXPEQ, NOUT, chease.log, ...)
- run_dir_top/             small files from the scratch run directory
- errors/                  error digest, console.log, logs of other attempts

Scratch run directory: {run_dir}
"""
    (stage_dir / "README.md").write_text(txt, encoding="utf-8")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--case-dir", type=Path, default=DEFAULT_CASE)
    ap.add_argument("--run-dir", type=Path, default=None,
                    help="scratch run directory (default: run_dir in cheasebs_acceptance.json)")
    ap.add_argument("--source-gfile", type=Path, default=None,
                    help="source EFIT (default: g*.00700 in the case dir)")
    ap.add_argument("--baseline-gfile", type=Path, default=DEFAULT_BASELINE,
                    help="unscaled cheaseBS EQDSK for comparison (test K Te_ped_scale_1.000)")
    ap.add_argument("--out", type=Path, default=None,
                    help="where the tarball goes (default: $SCRATCH, else the test directory's parent)")
    ap.add_argument("--name", default=None, help="tarball base name")
    ap.add_argument("--analysis-only", action="store_true", help="write review/ and stop")
    ap.add_argument("--skip-ogyropsi", action="store_true", help="drop ogyropsi.dat files (~24 MB each)")
    ap.add_argument("--with-baseline", action="store_true", help="also pack the run's cheasebs_baseline/")
    args = ap.parse_args()

    case = args.case_dir.resolve()
    test_dir = case.parent.parent
    acc = load_json(case / "cheasebs_acceptance.json")
    conv = load_json(case / "convergence_summary.json")
    run_dir = args.run_dir or (Path(acc["run_dir"]) if acc.get("run_dir") else None)
    source_g = args.source_gfile or next(
        (p for p in sorted(case.glob("g*")) if not p.name.endswith(".geqdsk")), None)
    if source_g is None:
        sys.exit(f"no source EFIT g-file in {case}")

    review_dir = test_dir / "review"
    C = run_analysis(case, source_g, args.baseline_gfile, run_dir, review_dir)
    tb = error_digest(case, test_dir, run_dir, review_dir / "errors")
    counts = C.counts()
    print(f"checks: {counts}  ->  {review_dir / 'check_report.md'}")
    for r in C.rows:
        if r["status"] in ("WARN", "FAIL"):
            print(f"  {r['status']:4s} {r['check']}: {_fmt(r['value'])} (expected {_fmt(r['expected'])}) {r['note']}")
    if tb:
        print(f"console.log: {len(tb)} failure line(s) from earlier attempts, see errors/error_digest.md")
    if args.analysis_only:
        return

    import os
    out = (args.out or Path(os.environ.get("SCRATCH", test_dir.parent))).resolve()
    name = args.name or f"{case.parent.name}_{test_dir.name}_review_{datetime.date.today():%Y%m%d}"
    stage_dir = out / name
    if stage_dir.exists():
        shutil.rmtree(stage_dir)
    copied = stage(case, test_dir, run_dir, review_dir, stage_dir, args.skip_ogyropsi, args.with_baseline, conv)
    shutil.move(str(stage_dir / "review" / "errors"), str(stage_dir / "errors"))
    write_readme(stage_dir, copied, run_dir, counts)

    tar_path = out / f"{name}.tar.gz"
    with tarfile.open(tar_path, "w:gz") as tf:
        tf.add(stage_dir, arcname=name)
    shutil.rmtree(stage_dir)
    print(f"tarball: {tar_path}  ({tar_path.stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
