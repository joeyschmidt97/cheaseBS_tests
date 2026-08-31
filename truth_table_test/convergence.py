"""The convergence path, not just the endpoint.

Scoring one iteration answers "is this equilibrium sound". It does not answer
"did the outer loop get there, or stop somewhere on the way", and on the NSTX
`I*` branch those are different questions: cheaseBS's own documentation records
iteration 00 as the validated replay point and later outer iterations as
numerically stable but drifting away from it. A run can therefore be perfectly
physical at iteration 00, perfectly physical at iteration 12, and disagree with
itself between the two.

So three things get read here, off `iteration_log.csv` and the per-iteration
EQDSKs:

* **residuals against their own tolerances** -- `tol_ip_rel`, `tol_bs`, `tol_q`,
  `tol_a` as the run was configured, not as this checker guesses;
* **the shape of the trace** -- falling, saturated, or oscillating. A loop can
  saturate at a value that never meets its tolerance, which is what cheaseBS
  does with `tol_ip_rel` under `qspec=on`, and that is a different failure from
  a residual that is still ringing;
* **iteration 00 versus the final iteration** -- the drift the NSTX branch is
  documented to have, measured rather than assumed, so a result can say which
  iteration it is quoting and how far apart the two are.
"""

from __future__ import annotations

import csv as _csv
import os

import numpy as np

from truth_table import (FAIL, FLAG, PASS, SKIP, Equilibrium, _grade, _rel,
                         _iteration_eqdsk, row)

# Residual column -> the config key that judges it, and a default for runs whose
# config did not record one. Defaults are the campaign settings in
# NSTXU_lithium_study/reshape_convergence/reshape_helpers.py (CHEASEBS).
RESIDUALS = (
    ("ip_error_rel", "tol_ip_rel", 0.02, "Ip closure"),
    ("bootstrap_change_rel", "tol_bs", 1e-3, "bootstrap change between iterations"),
    ("q_change_rel", "tol_q", 1e-3, "q change between iterations"),
    ("amplitude_change_rel", "tol_a", 1e-4, "driven-amplitude change"),
)


def read_iteration_log(path):
    """iteration_log.csv into {column: array}, non-numeric cells as NaN.

    Column names come from the cheaseBS driver, so nothing is assumed about the
    schema beyond the names being stable; every column that parses numerically
    is kept and the caller selects by name.
    """
    if not path or not os.path.isfile(path):
        return {}
    with open(path, newline="") as fh:
        rows = list(_csv.DictReader(fh))
    if not rows:
        return {}
    out = {}
    for key in rows[0]:
        vals = []
        for r in rows:
            try:
                vals.append(float((r.get(key) or "").strip()))
            except (TypeError, ValueError):
                vals.append(np.nan)
        arr = np.asarray(vals, dtype=float)
        if np.isfinite(arr).any():
            out[key] = arr
    return out


def _finite(arr):
    return arr[np.isfinite(arr)] if arr is not None else np.asarray([])


def saturation_iteration(values, rel_tol=0.05):
    """First iteration past which a metric stops moving by more than rel_tol.

    Saturation is the practically useful number -- the point after which more
    iterations buy nothing -- and is distinct from convergence, which is a
    comparison against a tolerance.
    """
    v = _finite(np.asarray(values, dtype=float))
    if v.size < 3:
        return None
    final = v[-1]
    scale = max(abs(final), 1e-30)
    for i in range(v.size):
        if np.all(np.abs(v[i:] - final) / scale <= rel_tol):
            return int(i)
    return None


def _trend(v):
    """(direction, sign_flips) for the tail of a residual trace.

    direction is the ratio last/first over the finite tail: below 1 means the
    residual fell. sign_flips counts direction changes in the successive
    differences of the last six points -- ringing, as distinct from a monotone
    approach or a flat saturation.
    """
    v = _finite(np.asarray(v, dtype=float))
    if v.size < 2:
        return np.nan, 0
    ratio = abs(v[-1]) / max(abs(v[0]), 1e-30)
    tail = v[-min(6, v.size):]
    d = np.diff(tail)
    d = d[np.abs(d) > 1e-12 * max(np.abs(tail).max(), 1e-30)]
    flips = int(np.sum(np.sign(d[1:]) != np.sign(d[:-1]))) if d.size > 1 else 0
    return float(ratio), flips


def build_rows(files, T):
    """Convergence-path rows for one run directory.

    T is the resolved tolerance dict from checks.build_rows, so the Ip bound is
    the same number in both blocks.
    """
    cfg, summ = files["cfg"], files["summary"]
    log = read_iteration_log(files.get("iteration_log"))
    rows = []
    add = rows.append

    if not log:
        add(row("iteration_trace", "convergence", "present", "iteration_log.csv",
                "absent", verdict=SKIP, fmt="{}",
                note="no iteration_log.csv: the endpoint can be scored but the "
                     "path to it cannot"))
        return rows

    n_iter = int(_finite(log.get("iteration", np.asarray([]))).size)
    max_iter = cfg.get("max_iter")
    at_cap = bool(max_iter and n_iter >= int(max_iter))
    add(row("iterations_run", "convergence", "not truncated",
            "max_iter " + str(max_iter), n_iter, None, None,
            FLAG if at_cap else PASS, fmt="{}",
            note="hit the cap: a residual that is still falling was cut off, so "
                 "'slow' and 'diverging' cannot be told apart here"
                 if at_cap else ""))

    for col, tol_key, tol_default, label in RESIDUALS:
        v = log.get(col)
        if v is None or not np.isfinite(v).any():
            add(row(col.replace("_rel", ""), "convergence", "<= tol", None, None,
                    verdict=SKIP,
                    note="column absent or all-NaN (normal for the first "
                         "iteration of a lagged quantity)"))
            continue
        tol = float(cfg.get(tol_key, tol_default))
        # Ip is judged on the shared tolerance so the two blocks cannot disagree.
        if tol_key == "tol_ip_rel":
            tol = T["ip_rel"]
        fin = _finite(v)
        final = float(fin[-1])
        ratio, flips = _trend(v)
        sat = saturation_iteration(v)
        verdict = _grade(final, tol, 5 * tol)
        note = "%s. first %.3g -> final %.3g" % (label, fin[0], final)
        if sat is not None:
            note += ", saturated at iteration %d" % sat
        if flips >= 3:
            note += ("; %d sign flips in the tail -- ringing, not converging. "
                     "Try a smaller bootstrap_mix / istar_mix before reading "
                     "physics off the final iteration." % flips)
            verdict = FLAG if verdict == PASS else verdict
        elif np.isfinite(ratio) and ratio > 1.0 and verdict != PASS:
            note += "; residual grew from its first value -- diverging"
            verdict = FAIL
        add(row(col.replace("_rel", ""), "convergence", "<= tol", 0.0, final,
                final, tol, verdict, note=note))

    # The workflow's own verdict, kept separate from the residual rows: the
    # acceptance policy deliberately does not gate on it.
    conv = summ.get("final_converged")
    add(row("driver_converged_flag", "convergence", "report only",
            "tol set above", str(conv), None, None,
            PASS if conv else FLAG if conv is False else SKIP, fmt="{}",
            note="cheaseBS's own flag. False with all residual rows passing means "
                 "one unmet tolerance, usually tol_ip_rel under qspec=on; False "
                 "with a FLAGged residual is the real thing"))

    # Iteration 00 versus the final iteration: the documented NSTX drift,
    # measured. Geometry is skipped -- header Ip, axis pressure and edge q
    # answer this and the contour integrals are the expensive part.
    its = files.get("iterations") or []
    if len(its) >= 2:
        e0, ef = _iteration_eqdsk(its[0][1]), _iteration_eqdsk(its[-1][1])
        if e0 and ef and os.path.realpath(e0) != os.path.realpath(ef):
            q0 = Equilibrium(e0, want_geometry=False)
            qf = Equilibrium(ef, want_geometry=False)
            drift = {
                "Ip": _rel(qf.ip, q0.ip),
                "p_axis": _rel(qf.p[0], q0.p[0]),
                "q(0.985)": _rel(qf.q_at(0.985, "rho_pol"), q0.q_at(0.985, "rho_pol")),
            }
            worst = max(abs(v) for v in drift.values() if np.isfinite(v))
            add(row("iter00_vs_final_drift", "convergence", "small", 0.0, worst,
                    worst, 0.01, _grade(worst, 0.01, 0.10),
                    note="max relative difference between iteration %d and "
                         "iteration %d (%s). cheaseBS documents iteration 00 as "
                         "the validated NSTX I* replay point and later outer "
                         "iterations as drifting, so quote one explicitly; a "
                         "large value here is the drift, not a rejection."
                         % (its[0][0], its[-1][0],
                            ", ".join("%s %+.2f%%" % (k, 100 * v)
                                      for k, v in drift.items()
                                      if np.isfinite(v)))))
        else:
            add(row("iter00_vs_final_drift", "convergence", "small", None, None,
                    verdict=SKIP,
                    note="only one distinct per-iteration EQDSK on disk"))
    else:
        add(row("iter00_vs_final_drift", "convergence", "small", None, None,
                verdict=SKIP,
                note="fewer than two iteration directories: nothing to compare"))
    return rows
