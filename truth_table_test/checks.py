"""The truth table itself: one row per quantity, PASS / FLAG / FAIL / SKIP.

Two blocks, and they fail in opposite directions.

**Controls** are imposed by the run and must survive the profile scan: the
boundary, the toroidal field, the target current, the enclosed volume, the edge
shaping, and the reference fast-pressure / driven-current decomposition. A
change here means the comparison is not a profile scan at all.

**Responses** are what the scaled pfiles are supposed to drive: thermal
pressure, axis pressure, beta, stored energy, the magnetic axis and internal
flux-surface centre, the bootstrap current and the replayed total I*. A
quantity in this block that does *not* move is the failure -- specifically the
signature of profiles that never reached CHEASE.

q sits between the two and is deliberately not gated on fidelity: a reshaped
pedestal drives more bootstrap current and is supposed to move q. Only the
blowup bound is enforced. See CheasebsAcceptance in TPED's cheasebs_runner for
the same split.

Quantity definitions and readers live in truth_table.py.
"""

from __future__ import annotations

import os

import numpy as np

import convergence

from truth_table import (FAIL, FLAG, PASS, SKIP, SPECIES, Equilibrium,
                         _direction, _grade, _rel, _trapz,
                         boundary_displacement, read_csv_columns, read_xy, row,
                         thermal_pressure)

DEFAULT_TOL = dict(
    ip_rel=0.02,                            # overridden by the run's own tol_ip_rel
    bt_rel=1e-6,
    lcfs_warn=0.005, lcfs_fail=0.02,        # fraction of a_minor
    shape_warn=0.02, shape_fail=0.05,       # kappa_95 relative, delta_95 absolute
    vol_warn=0.02, vol_fail=0.05,
    axis_noise_m=1e-4,                      # 0.1 mm: below this it is grid noise
    q_edge_doc=0.01, q_edge_fail=0.30,
    q_blowup=0.5,
    istar_closure_warn=1e-9, istar_closure_fail=1e-6,
    p_axis_resid_warn=0.05, p_axis_resid_fail=0.15,
    dead=1e-3,                              # the "it never moved" floor
)


def _pth_reference(files):
    """Reference thermal pressure: the REF_ profiles, or the baseline CSV.

    The REF_ profiles are the honest reference -- they are the arrays cheaseBS
    was handed. The baseline decomposition's p_th_pa column is the same physics
    one interpolation later, and is the only reference available for runs that
    predate reference_profiles/, so it is used as a labelled fallback.
    """
    rho, p, used = thermal_pressure(files["profiles_before"])
    if rho is not None:
        return rho, p, "REF_profiles (" + "".join(used) + ")"
    if files.get("baseline_csv"):
        col = read_csv_columns(files["baseline_csv"])
        if "rhot" in col and "p_th_pa" in col:
            return col["rhot"], col["p_th_pa"], "baseline profiles.csv p_th_pa"
    return None, None, None


def ratio_stats(rho_a, a, rho_b, b, radii=None, rho_max=0.995):
    """Where and by how much a/b changed, on the reference grid.

    A pedestal-height scale changes a narrow radial band, so a whole-profile
    median understates it and can even invert: measured on the 2026-08-28
    campaign, 132543 Te at scale 0.70 gives a median ratio of 1.013 while the
    ratio at rho_tor 0.9 is 0.927 -- the median is dominated by the core, which
    moves slightly the other way when the Stefanikova fit renormalises. Scoring
    the median therefore failed eight Te points that had scaled correctly.

    So the scored number is the ratio at the radius the scan is *about*: the
    GENE analysis radius when the campaign recorded one, else rho_tor 0.9. The
    median and the peak deviation come back alongside it as context, not as the
    verdict.

    Returns {at_radius, radius, median, peak, peak_rho}.
    """
    out = dict(at_radius=np.nan, radius=np.nan, median=np.nan,
               peak=np.nan, peak_rho=np.nan)
    m = (rho_b <= rho_max) & np.isfinite(b) & (np.abs(b) > 0)
    if m.sum() < 5:
        return out
    o = np.argsort(rho_a)
    rho = rho_b[m]
    ratio = np.interp(rho, rho_a[o], a[o]) / b[m]
    ok = np.isfinite(ratio)
    if ok.sum() < 5:
        return out
    rho, ratio = rho[ok], ratio[ok]

    # The scoring radius: the outermost analysis radius inside the fitted range,
    # since that is the one nearest the pedestal the scan moved.
    cand = [float(r) for r in (radii or []) if np.isfinite(float(r))]
    cand = [r for r in cand if rho.min() <= r <= rho.max()]
    radius = max(cand) if cand else float(min(0.9, rho.max()))

    i = int(np.nanargmax(np.abs(ratio - 1.0)))
    out.update(at_radius=float(np.interp(radius, rho, ratio)),
               radius=radius,
               median=float(np.nanmedian(ratio)),
               peak=float(ratio[i]), peak_rho=float(rho[i]))
    return out


def build_rows(files, tol=None):
    """Every truth-table row for one run directory."""
    T = dict(DEFAULT_TOL)
    if tol:
        T.update(tol)
    cfg, summ, acc = files["cfg"], files["summary"], files["acceptance_rec"]
    scale = files["case"]["scale"]
    if cfg.get("tol_ip_rel"):
        T["ip_rel"] = float(cfg["tol_ip_rel"])

    rows = []
    add = rows.append

    # ================= controls: must not move =================

    if files["eqdsk_after"] is None:
        add(row("run_completed", "control", "reconstruction exists", "EQDSK*.OUT",
                "missing", verdict=FAIL, fmt="{}",
                note="no EQDSK*.OUT: cheaseBS produced no equilibrium here"))
        return rows
    add(row("run_completed", "control", "reconstruction exists", "EQDSK*.OUT",
            os.path.basename(files["eqdsk_after"]), verdict=PASS, fmt="{}",
            note="iteration " + str(files["iteration_index"])))

    # Solver state. Not converging is not automatically wrong -- the acceptance
    # policy deliberately does not gate on cheaseBS's own flag, which is
    # permanently False on 132588's usable equilibrium -- but a run that stopped
    # by exhausting max_iter has to be visible, or "slow" and "diverging" cannot
    # be told apart.
    iters, conv = summ.get("iterations"), summ.get("final_converged")
    max_iter = cfg.get("max_iter")
    at_cap = bool(iters and max_iter and iters >= int(max_iter))
    add(row("solver_state", "control", "completed, not truncated",
            "max_iter " + str(max_iter),
            str(iters) + " iters, converged=" + str(conv), fmt="{}",
            verdict=(FLAG if (at_cap or conv is False) else PASS if conv else SKIP),
            note=("hit the iteration cap" if at_cap else
                  "cheaseBS converged flag is False; the acceptance policy does not "
                  "gate on it, but read the iteration trace" if conv is False else
                  "no convergence_summary.json" if not summ else "")))

    eq_a = Equilibrium(files["eqdsk_after"])
    eq_b = Equilibrium(files["eqdsk_before"]) if files["eqdsk_before"] else None
    if eq_b is None:
        add(row("source_eqdsk", "control", "present", "g<shot>.<time>", "missing",
                verdict=FAIL, fmt="{}",
                note="no source EQDSK in the run dir: nothing to compare against. "
                     "Pass --source-eqdsk."))

    target_ip = (summ.get("target_ip_a") or cfg.get("target_ip_a")
                 or (eq_b.ip if eq_b else None))
    ip_err = _rel(eq_a.ip, target_ip)
    add(row("ip_closure", "control", "at target", target_ip, eq_a.ip, ip_err,
            T["ip_rel"], _grade(ip_err, T["ip_rel"], 2 * T["ip_rel"]),
            note="Ip from the reconstruction's own header"))

    if eq_b:
        d = _rel(eq_a.b0, eq_b.b0)
        add(row("bt_invariant", "control", "unchanged", eq_b.b0, eq_a.b0, d,
                T["bt_rel"], _grade(d, T["bt_rel"], 1e-3),
                note="|Bctr| from the EQDSK headers"))

        # The fixed-boundary control, and the headline one.
        disp = boundary_displacement(eq_b, eq_a)
        a_min = eq_b.a_minor
        if disp and np.isfinite(a_min) and a_min > 0:
            rms, mx = disp["rms_m"] / a_min, disp["max_m"] / a_min
            add(row("lcfs_rms/a", "control", "unchanged", 0.0, rms, rms,
                    T["lcfs_warn"], _grade(rms, T["lcfs_warn"], T["lcfs_fail"]),
                    note="rms %.2f mm, a = %.3f m; sampled at common poloidal "
                         "angle about the source axis" % (disp["rms_m"] * 1e3, a_min)))
            add(row("lcfs_max/a", "control", "unchanged (report)", 0.0, mx, mx,
                    T["lcfs_fail"],
                    PASS if mx <= T["lcfs_fail"] else FLAG,
                    note="max %.2f mm. The X-point region dominates this number "
                         "and a hard bound there would reject sound "
                         "reconstructions, so it flags rather than fails; "
                         "lcfs_rms/a is the gate." % (disp["max_m"] * 1e3)))
        else:
            add(row("lcfs_rms/a", "control", "unchanged", None, None, verdict=SKIP,
                    note="one of the EQDSKs carries no boundary polygon"))

        sa, sb = eq_a.shape_at(0.95), eq_b.shape_at(0.95)
        if sa and sb:
            dk = _rel(sa["kappa"], sb["kappa"])
            add(row("kappa_95", "control", "unchanged", sb["kappa"], sa["kappa"],
                    dk, T["shape_warn"],
                    _grade(dk, T["shape_warn"], T["shape_fail"])))
            dd = sa.get("delta", np.nan) - sb.get("delta", np.nan)
            add(row("delta_95", "control", "unchanged", sb.get("delta"),
                    sa.get("delta"), dd, T["shape_warn"],
                    _grade(dd, T["shape_warn"], T["shape_fail"]),
                    note="absolute change, mean of upper and lower"))

        dv = _rel(eq_a.volume, eq_b.volume)
        add(row("plasma_volume", "control", "unchanged", eq_b.volume, eq_a.volume,
                dv, T["vol_warn"], _grade(dv, T["vol_warn"], T["vol_fail"]),
                note="m^3, Pappus integral over closed contours. The 0.02 bound "
                     "is set by the grid change alone: two identity replays "
                     "(132588, 129015) move the volume 0.8-1.2% with the "
                     "boundary held fixed, because CHEASE writes a 513x513 "
                     "EQDSK and the source EFIT is 129x129."))

    # The fast-pressure and driven-current decomposition must be the REFERENCE
    # one. Aliasing the reference profiles onto the scaled ones rebuilds the
    # baseline from the scaled profiles, which floors p_total at the scaled
    # pressure and discards every downward scan -- the 2026-08-26 bug.
    if files["profiles_aliased"]:
        add(row("pfast_reference_frozen", "control", "reference decomposition",
                "REF_profiles distinct", "aliased onto run profiles",
                verdict=FAIL, fmt="{}",
                note="reference profiles ARE the run profiles (" +
                     ",".join(files["profiles_aliased"]) + "): p_fast moved with "
                     "the scan, so this point is not a thermal-profile-only scan"))
    elif all(files["profiles_before"].get(s) for s, _ in SPECIES):
        note = ""
        if files.get("baseline_csv"):
            pf = read_csv_columns(files["baseline_csv"]).get("p_fast_pa")
            if pf is not None and np.isfinite(pf).any():
                note = ("baseline p_fast max %.1f Pa (compared across cases in "
                        "--campaign mode)" % np.nanmax(pf))
        add(row("pfast_reference_frozen", "control", "reference decomposition",
                "REF_profiles distinct", "ok", verdict=PASS, note=note, fmt="{}"))
    else:
        add(row("pfast_reference_frozen", "control", "reference decomposition",
                "REF_profiles", "absent", verdict=SKIP, fmt="{}",
                note="no reference_profiles/: predates the reference-profile fix, "
                     "or an identity replay"))

    if files.get("baseline_istar_csv"):
        col = read_csv_columns(files["baseline_istar_csv"])
        tot, rec = col.get("total_istar"), col.get("reconstructed_total_istar")
        if tot is not None and rec is not None:
            den = max(np.nanmax(np.abs(tot)), 1e-30)
            resid = float(np.nanmax(np.abs(rec - tot)) / den)
            add(row("istar_baseline_closure", "control", "machine precision", 0.0,
                    resid, resid, T["istar_closure_warn"],
                    _grade(resid, T["istar_closure_warn"], T["istar_closure_fail"]),
                    note="max|reconstructed - total| / max|total| on the baseline grid"))

    # ================= responses: must move, and the right way =================

    rho_ref, pth_ref, ref_label = _pth_reference(files)
    rho_act, pth_act, used_act = thermal_pressure(files["profiles_after"])
    radii = files["case"].get("radii")

    stats = None
    scan_active = True
    if rho_ref is not None and rho_act is not None:
        stats = ratio_stats(rho_act, pth_act, rho_ref, pth_ref, radii=radii)
        # Gate the rest of the response block on the scored radius, not the
        # median: a real pedestal scan can leave the median at 1.0.
        scan_active = bool(np.isfinite(stats["at_radius"])
                           and abs(stats["at_radius"] - 1.0) > T["dead"])
    elif scale in (None, 1.0):
        # No measurable profile change and no scan direction: treat the run as
        # an identity replay and report the responses instead of scoring them.
        scan_active = False

    if stats is not None and np.isfinite(stats["at_radius"]):
        v, why = _direction(stats["at_radius"] - 1.0, scale, T["dead"],
                            gated=not (scale in (None, 1.0) and not scan_active))
        add(row("p_th_at_radius", "response", "changes", 1.0,
                stats["at_radius"], stats["at_radius"] - 1.0, None, v,
                note=("p_th ratio at rho_tor %.3f (%s). whole-profile median "
                      "%.3f, peak %.3f at rho_tor %.3f -- the median understates "
                      "a pedestal-only scale and can invert, so it is context "
                      "only. reference = %s; active species %s. %s"
                      % (stats["radius"],
                         "campaign analysis radius" if radii else "default 0.9",
                         stats["median"], stats["peak"], stats["peak_rho"],
                         ref_label, "".join(used_act), why)).strip()))
    else:
        add(row("p_th_at_radius", "response", "changes", None, None,
                verdict=SKIP, note="reference or active profiles unavailable"))

    # The direct pfile-ingestion test: p_total = p_th(active) + p_fast(reference).
    if files.get("baseline_csv") and rho_act is not None:
        col = read_csv_columns(files["baseline_csv"])
        pf, rhot = col.get("p_fast_pa"), col.get("rhot")
        if pf is not None and rhot is not None:
            o = np.argsort(rhot)
            pred = float(np.interp(0.0, rhot[o], pf[o])) + float(
                np.interp(0.0, rho_act, pth_act))
            got = float(eq_a.p[0])
            resid = _rel(got, pred)
            add(row("p_axis_ingestion", "response", "matches inputs", pred, got,
                    resid, T["p_axis_resid_warn"],
                    _grade(resid, T["p_axis_resid_warn"], T["p_axis_resid_fail"]),
                    note="axis p of the reconstruction vs p_fast(reference) + "
                         "p_th(active), both in Pa"))

    if eq_b:
        d = _rel(eq_a.p[0], eq_b.p[0])
        v, why = _direction(d, scale, T["dead"], gated=scan_active)
        add(row("p_axis_response", "response", "changes", eq_b.p[0], eq_a.p[0],
                d, None, v, note=("Pa, on axis. " + why).strip()))

        d = _rel(eq_a.beta_t, eq_b.beta_t)
        v, why = _direction(d, scale, T["dead"], gated=scan_active)
        add(row("beta_t", "response", "changes", eq_b.beta_t, eq_a.beta_t, d,
                None, v, note=("2 mu0 <p>_V / B0^2. " + why).strip()))
        add(row("beta_N", "response", "report only", eq_b.beta_n, eq_a.beta_n,
                _rel(eq_a.beta_n, eq_b.beta_n), None, PASS,
                note="beta_t[%] a B0 / Ip[MA]; reported, not gated"))

        d = _rel(eq_a.w_thermal, eq_b.w_thermal)
        v, why = _direction(d, scale, T["dead"], gated=scan_active)
        add(row("W_stored", "response", "changes", eq_b.w_thermal, eq_a.w_thermal,
                d, None, v,
                note=("J, 3/2 int p dV on the EQDSK's own total p -- thermal plus "
                      "the fixed fast-ion part. " + why).strip()))

        # Internal geometry: the reshape a fixed boundary hides.
        dR = eq_a.r_axis - eq_b.r_axis
        moved = np.isfinite(dR) and abs(dR) > T["axis_noise_m"]
        add(row("axis_R_shift_mm", "response", "changes", eq_b.r_axis * 1e3,
                eq_a.r_axis * 1e3, dR * 1e3, T["axis_noise_m"] * 1e3,
                SKIP if not np.isfinite(dR) else PASS if moved else FLAG,
                note="" if moved else "below the 0.1 mm grid-noise floor: not a "
                                      "resolved response"))
        add(row("axis_Z_shift_mm", "response", "report only", eq_b.z_axis * 1e3,
                eq_a.z_axis * 1e3, (eq_a.z_axis - eq_b.z_axis) * 1e3, None, PASS,
                note="vertical axis motion at fixed boundary; large values are suspect"))

        ra, rb = eq_a.r_geo_at(0.5), eq_b.r_geo_at(0.5)
        d = (ra - rb) * 1e3 if np.isfinite(ra) and np.isfinite(rb) else np.nan
        add(row("shafranov_rho0.5_mm", "response", "changes",
                rb * 1e3 if np.isfinite(rb) else None,
                ra * 1e3 if np.isfinite(ra) else None, d,
                T["axis_noise_m"] * 1e3,
                SKIP if not np.isfinite(d) else
                PASS if abs(d) > T["axis_noise_m"] * 1e3 else FLAG,
                note="flux-surface centre at rho_tor 0.5; outward at higher pressure"))

    q_err = acc.get("q_errors_rel") or {}
    if q_err:
        worst = max(abs(v) for v in q_err.values())
        add(row("q_at_analysis_radii", "response", "may change, bounded",
                "source q", worst, worst, T["q_blowup"],
                _grade(worst, T["q_blowup"], T["q_blowup"]),
                note="max |dq|/q at the GENE analysis radii ("
                     + ", ".join("%s: %+.2f%%" % (k, 100 * v)
                                 for k, v in q_err.items())
                     + "). Gated only as a blowup bound -- a reshaped pedestal is "
                       "supposed to move q."))

    q_edge = acc.get("q_edge_error_rel")
    if q_edge is None:
        q_edge = summ.get("q_report_rel_diff")
    if q_edge is not None:
        add(row("q_edge_error", "response", "near reference", 0.0, q_edge, q_edge,
                T["q_edge_doc"], _grade(q_edge, T["q_edge_doc"], T["q_edge_fail"]),
                note="1% is the documented QSPEC figure; these NSTX cases carry a "
                     "standing 10-20% offset on the outermost surface (grid/X-point "
                     "artifact), so >1% is a FLAG to read, not a rejection"))

    if eq_b:
        add(row("q0_qmin_q95", "response", "report only",
                "%.3g/%.3g/%.3g" % (eq_b.q[0], np.abs(eq_b.q).min(),
                                    eq_b.q_at(0.95, "rho_pol")),
                "%.3g/%.3g/%.3g" % (eq_a.q[0], np.abs(eq_a.q).min(),
                                    eq_a.q_at(0.95, "rho_pol")),
                None, None, PASS, fmt="{}",
                note="q0 / |q|min / q95(rho_pol); read for discontinuities or a "
                     "rational-surface crossing"))

    # Bootstrap current must respond to the changed gradients and collisionality.
    bs_case = None
    if files["iteration_dir"]:
        p = os.path.join(files["iteration_dir"], "bootstrap_diagnostic_rhop.txt")
        bs_case = p if os.path.isfile(p) else None
    if bs_case and files.get("baseline_csv"):
        col = read_csv_columns(files["baseline_csv"])
        rho_r = col.get("rhop")
        j_r = col.get("bootstrap_parallel_current_density")
        rho_c, j_c = read_xy(bs_case)
        if rho_r is not None and j_r is not None:
            o, oc = np.argsort(rho_r), np.argsort(rho_c)
            i_ref = float(_trapz(np.abs(j_r[o]), rho_r[o]))
            i_case = float(_trapz(np.abs(j_c[oc]), rho_c[oc]))
            d = _rel(i_case, i_ref)
            v, why = _direction(d, scale, T["dead"], gated=scan_active)
            add(row("bootstrap_response", "response", "changes", i_ref, i_case,
                    d, None, v,
                    note=("int |j_bs| drhop: reference decomposition vs this "
                          "iteration's diagnostic. A shape proxy, not I_bs in A. "
                          + why).strip()))
    else:
        add(row("bootstrap_response", "response", "changes", None, None,
                verdict=SKIP,
                note="no bootstrap_diagnostic_rhop.txt or no baseline profiles.csv"))

    # The replayed total I* is what CHEASE actually solved with.
    istar_case = None
    if files["iteration_dir"]:
        p = os.path.join(files["iteration_dir"], "istar_target_rhop.txt")
        istar_case = p if os.path.isfile(p) else None
    if istar_case and files.get("baseline_istar_csv"):
        col = read_csv_columns(files["baseline_istar_csv"])
        rho_r, i_r = col.get("rho"), col.get("total_istar")
        rho_c, i_c = read_xy(istar_case)
        if rho_r is not None and i_r is not None:
            o, oc = np.argsort(rho_r), np.argsort(rho_c)
            i_on_ref = np.interp(rho_r[o], rho_c[oc], i_c[oc])
            den = max(np.nanmax(np.abs(i_r)), 1e-30)
            rms = float(np.sqrt(np.nanmean((i_on_ref - i_r[o]) ** 2)) / den)
            add(row("istar_total_response", "response", "changes", 0.0, rms, rms,
                    T["dead"],
                    SKIP if not np.isfinite(rms) else
                    PASS if rms > T["dead"] else FAIL,
                    note="rms(I*_target - I*_baseline)/max|I*_baseline|. Zero means "
                         "the replayed current never saw the profile change."))
    else:
        add(row("istar_total_response", "response", "changes", None, None,
                verdict=SKIP, note="no istar_target_rhop.txt or no baseline I* CSV"))

    bad = []
    if not np.isfinite(eq_a.q).all():
        bad.append("q has non-finite values")
    elif np.nanmin(np.abs(eq_a.q)) <= 0:
        bad.append("|q| reaches zero")
    if not np.isfinite(eq_a.p).all():
        bad.append("p has non-finite values")
    elif np.nanmin(eq_a.p) < -1.0:
        bad.append("p goes negative (%.3g Pa)" % np.nanmin(eq_a.p))
    if eq_a.geom is None:
        bad.append("flux surfaces not contourable: %s" % eq_a.geom_error)
    add(row("physicality", "control", "physical", "q>0, p>=0, nested",
            "; ".join(bad) if bad else "ok", verdict=FAIL if bad else PASS,
            fmt="{}"))

    # How the loop got here, judged against the run's own tolerances.
    rows.extend(convergence.build_rows(files, T))

    return rows


# --------------------------------------------------------------------------
# cross-case checks: only a campaign can answer these
# --------------------------------------------------------------------------

def _value(rows, name):
    for r in rows:
        if r["check"] == name:
            return r
    return None


def campaign_rows(cases):
    """Checks that need more than one point: invariance across the scan, and
    low-below-high ordering per axis.

    cases : list of dicts with 'label', 'axis', 'scale', 'rows'.
    """
    out = []
    add = out.append

    # Every control that is supposed to be a run *setting* must be identical
    # across the scan, or the cases are not comparable to each other -- however
    # well each one scores on its own.
    for name, tol_rel in (("bt_invariant", 1e-9), ("ip_closure", None)):
        vals = [(c["label"], _value(c["rows"], name)) for c in cases]
        refs = [(lab, r["reference"]) for lab, r in vals
                if r and isinstance(r["reference"], (int, float))]
        if len(refs) < 2:
            continue
        base = refs[0][1]
        worst = max(abs((v - base) / base) if base else np.inf for _, v in refs)
        label = {"bt_invariant": "Bt target identical",
                 "ip_closure": "Ip target identical"}[name]
        add(row(label, "campaign", "identical", base, refs[-1][1], worst,
                tol_rel or 1e-9,
                PASS if worst <= (tol_rel or 1e-9) else FAIL,
                note="a target that moves between cases invalidates the "
                     "comparison, whatever each case scores alone"))

    # Ordering: for each scan axis, the low point must sit below the high point
    # in every pressure-driven quantity. A single case can look plausible while
    # the pair is ordered backwards.
    by_axis = {}
    for c in cases:
        if c["axis"] is None or c["scale"] is None:
            continue
        by_axis.setdefault(c["axis"], []).append(c)
    for axis, group in sorted(by_axis.items()):
        group = sorted(group, key=lambda c: c["scale"])
        if len(group) < 2:
            continue
        lo, hi = group[0], group[-1]
        for name in ("p_th_at_radius", "p_axis_response", "beta_t",
                     "W_stored", "bootstrap_response"):
            rlo, rhi = _value(lo["rows"], name), _value(hi["rows"], name)
            if not (rlo and rhi):
                continue
            vlo, vhi = rlo["value"], rhi["value"]
            if not (isinstance(vlo, (int, float)) and isinstance(vhi, (int, float))
                    and np.isfinite(vlo) and np.isfinite(vhi)):
                continue
            ok = vlo < vhi
            add(row("%s/%s" % (axis.replace("_ped_scale", ""), name), "campaign",
                    "low < high", "%.4g (x%.2f)" % (vlo, lo["scale"]),
                    "%.4g (x%.2f)" % (vhi, hi["scale"]), None, None,
                    PASS if ok else FAIL, fmt="{}",
                    note="" if ok else "ordering reversed: the scale was not "
                                       "applied as intended, or another setting "
                                       "changed between the two points"))
    return out


# --------------------------------------------------------------------------
# rendering
# --------------------------------------------------------------------------

_W = (30, 21, 15, 15, 11, 9, 8)
_HEAD = ("check", "expect", "reference", "case", "delta", "tol", "verdict")


def _cell(v, fmt):
    if v is None:
        return "--"
    if isinstance(v, str):
        return v
    if isinstance(v, (int, float)):
        if not np.isfinite(v):
            return "--"
        try:
            return fmt.format(v)
        except (TypeError, ValueError):
            return str(v)
    return str(v)


_GROUP_ORDER = ("control", "response", "convergence", "campaign")


def render(rows, title=None, notes=True, width=112):
    """The table as plain text: fixed columns, notes indented under the row.

    Rows are shown grouped rather than in construction order -- the control /
    response / convergence split is the point of the table, and physicality is
    built last but belongs with the controls.
    """
    rows = sorted(rows, key=lambda r: _GROUP_ORDER.index(r["group"])
                  if r["group"] in _GROUP_ORDER else len(_GROUP_ORDER))
    lines = []
    if title:
        lines.append(title)
    lines.append("  " + "".join(h.ljust(w) for h, w in zip(_HEAD, _W)))
    lines.append("  " + "-" * sum(_W))
    group = None
    for r in rows:
        if r["group"] != group:
            group = r["group"]
            lines.append("  [%s]" % group)
        cells = [r["check"],
                 r["expect"],
                 _cell(r["reference"], r["fmt"]),
                 _cell(r["value"], r["fmt"]),
                 _cell(r["delta"], "{:+.3g}"),
                 _cell(r["tol"], "{:.3g}"),
                 r["verdict"]]
        lines.append("  " + "".join(str(c)[:w - 1].ljust(w)
                                    for c, w in zip(cells, _W)))
        if notes and r["note"]:
            for chunk in _wrap(r["note"], width - 8):
                lines.append("        " + chunk)
    return "\n".join(lines)


def _wrap(text, width):
    words, line, out = text.split(), "", []
    for w in words:
        if len(line) + len(w) + 1 > width:
            out.append(line)
            line = w
        else:
            line = (line + " " + w).strip()
    if line:
        out.append(line)
    return out


def tally(rows):
    counts = {PASS: 0, FLAG: 0, FAIL: 0, SKIP: 0}
    for r in rows:
        counts[r["verdict"]] = counts.get(r["verdict"], 0) + 1
    return counts
