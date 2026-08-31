"""Truth table for one cheaseBS reconstruction: what must hold, what must move.

A profile-scaling run is only interpretable if the *controls* survived it. The
boundary, the toroidal field and the target current are imposed; the thermal
pressure, the bootstrap current, the replayed I* and the internal flux geometry
are supposed to respond. Those two groups fail in opposite directions, so a
single "did it converge" flag cannot separate a working reshape from a run that
never ingested the scaled profiles at all -- the 2026-08-26 reference-profile
bug looked like a converged run for weeks.

This module holds the readers and the derived quantities both groups are
scored on -- EQDSK equilibria, GENE profiles, the baseline decomposition, the
per-iteration current artifacts. The table itself is in checks.py and the CLI
in run_truth_table.py. Nothing here re-runs CHEASE; it reads what a completed
run left on disk.

WHAT COUNTS AS THE REFERENCE

Not a separate baseline solve. Every cheaseBS run directory already carries the
"before" side of its own comparison:

* the source EQDSK, copied in under its original ``g<shot>.<time>`` name;
* ``reference_profiles/REF_profiles_{e,i,z}`` -- the untransformed profiles the
  source EQDSK is consistent with;
* ``cheasebs_baseline/`` -- the fast-pressure and current decomposition built
  from those reference profiles, which must NOT move when the active profiles
  are scaled.

So a run is self-checking, and the checker can be pointed at one directory.

WHAT IS DELIBERATELY NOT GATED

q at the analysis radii is *expected* to move: a taller pedestal drives more
bootstrap current and changes q, which is the reason CHEASE is in the loop.
Only the blowup bound is gated. The edge-q error carries a known standing
offset (10-20% on these NSTX cases, a grid/X-point artifact) and is reported
against the 1% documentation figure as a FLAG, never a FAIL.

USAGE

    python run_truth_table.py --run  <savedir>
    python run_truth_table.py --campaign <runs/129015_20260824_20-28-57>

See README.md in this directory.
"""

from __future__ import annotations

import glob
import json
import os
import re
import sys

import numpy as np

MU0 = 4.0e-7 * np.pi

# numpy renamed trapz -> trapezoid in 2.0; NERSC still ships 1.x in places.
_trapz = getattr(np, "trapezoid", None) or np.trapz

# p[Pa] = n[1e19 m^-3] * T[keV] * this. Checked against the p_e_pa column of a
# baseline profiles.csv: 4.45182 * 0.810705 * 1602.176 = 5782.44 Pa, and the
# CSV says 5782.436.
PA_PER_1E19_KEV = 1.602176634e3

SPECIES = (("e", "electron"), ("i", "deuterium"), ("z", "carbon"))

PASS, FLAG, FAIL, SKIP = "PASS", "FLAG", "FAIL", "SKIP"


# --------------------------------------------------------------------------
# imports that live outside this repo
# --------------------------------------------------------------------------

def _bootstrap_tped():
    """Make TPED importable without assuming how it was installed.

    reshape_helpers.py imports it bare (installed / PYTHONPATH), while
    ip_consistency_check.py hardcodes the Windows checkout. Try the import
    first, then the usual checkout roots, and say plainly what to set if none
    of them work -- a stack trace three frames into xarray does not.
    """
    try:
        import TPED  # noqa: F401
        return
    except ImportError:
        pass
    cands = [os.environ.get("TPED_ROOT"),
             os.path.expanduser("~/git"),
             os.path.expanduser("~/projects"),
             "C:/Users/joesc/git",
             "/global/u1/j/joeschm/projects"]
    for root in cands:
        if root and os.path.isdir(os.path.join(root, "TPED")):
            sys.path.insert(0, root)
            return
    raise ImportError(
        "TPED not importable. Set TPED_ROOT to the directory *containing* the "
        "TPED checkout, e.g. TPED_ROOT=C:/Users/joesc/git")


def _gfile_reader():
    _bootstrap_tped()
    from TPED.projects.discharge_tools.src.filetypes.gfile_data import GFileData
    return GFileData


def _geometry():
    _bootstrap_tped()
    from TPED.projects.discharge_tools.src.equilibrium_geometry import flux_surface_geometry
    return flux_surface_geometry


# --------------------------------------------------------------------------
# locating what a run left on disk
# --------------------------------------------------------------------------

def _find_source_gfile(run_dir):
    """The frozen SOURCE geometry, not the reconstruction.

    The input EFIT is copied in under its own g<shot>.<time> name and cheaseBS
    writes its result as EQDSK*.OUT. Globbing g* loosely is the trap that makes
    every before/after delta come out zero.
    """
    cands = [p for p in glob.glob(os.path.join(run_dir, "g*"))
             if os.path.isfile(p) and re.match(r"^g\d+\.\d+", os.path.basename(p))]
    return sorted(cands)[0] if cands else None


def _iteration_dirs(run_dir):
    """iteration_NN directories, wherever the driver's output_dir landed.

    The TPED wrapper sets output_dir = savedir, so they sit directly in the run
    directory. A bare cheaseBS config can point output_dir at a `runs/`
    subdirectory instead (that is the layout under data/*_runs/ in this repo),
    so look one level down as well.
    """
    hits = sorted(glob.glob(os.path.join(run_dir, "iteration_*")))
    if not hits:
        hits = sorted(glob.glob(os.path.join(run_dir, "runs", "iteration_*")))
    out = []
    for p in hits:
        m = re.search(r"iteration_(\d+)$", p)
        if m and os.path.isdir(p):
            out.append((int(m.group(1)), p))
    return sorted(out)


def _iteration_eqdsk(itdir):
    """The EQDSK that iteration wrote, preferring the positive-sign COCOS copy."""
    for pat in ("artifacts/EQDSK_COCOS_02_POS*.OUT", "artifacts/EQDSK*.OUT",
                "EQDSK_COCOS_02_POS*.OUT", "EQDSK*.OUT"):
        hits = sorted(glob.glob(os.path.join(itdir, pat)))
        if hits:
            return hits[0]
    return None


def locate(run_dir, source_eqdsk=None, iteration="final"):
    """Everything the checks read, resolved once, with the misses recorded.

    iteration : "final" (default) or an integer / "00". The NSTX I* branch
        documents iteration 00 as its validated replay point and later outer
        iterations as still drifting, so both have to be scoreable.
    """
    run_dir = os.path.abspath(run_dir)
    f = {"run_dir": run_dir, "notes": []}

    cfg_path = os.path.join(run_dir, "cheasebs_run_config.json")
    cfg = {}
    if os.path.isfile(cfg_path):
        try:
            cfg = json.load(open(cfg_path))
            f["config"] = cfg_path
        except Exception as exc:
            f["notes"].append(f"cheasebs_run_config.json unreadable: {exc}")
    else:
        f["notes"].append("no cheasebs_run_config.json; paths inferred from the directory")
    f["cfg"] = cfg

    # Active and reference profiles. Config paths are absolute and point at the
    # machine the solve ran on, so they are used only when they still resolve.
    f["profiles_after"], f["profiles_before"] = {}, {}
    for spec, longname in SPECIES:
        after = cfg.get(f"{longname}_profile")
        before = cfg.get(f"reference_{longname}_profile")
        if not (after and os.path.isfile(after)):
            after = os.path.join(run_dir, f"profiles_{spec}")
        if not (before and os.path.isfile(before)):
            before = os.path.join(run_dir, "reference_profiles", f"REF_profiles_{spec}")
        f["profiles_after"][spec] = after if os.path.isfile(after) else None
        f["profiles_before"][spec] = before if os.path.isfile(before) else None

    # The bug this layout exists to catch: before 2026-08-26 the reference
    # profiles WERE the run profiles, the baseline decomposition was rebuilt
    # from the scaled profiles, and every downward scan was silently discarded.
    same = [s for s, _ in SPECIES
            if f["profiles_before"][s] and f["profiles_after"][s]
            and os.path.normcase(os.path.realpath(f["profiles_before"][s]))
            == os.path.normcase(os.path.realpath(f["profiles_after"][s]))]
    f["profiles_aliased"] = same

    src = source_eqdsk or cfg.get("eqdsk")
    if not (src and os.path.isfile(src)):
        src = _find_source_gfile(run_dir)
    f["eqdsk_before"] = src if src and os.path.isfile(src) else None

    its = _iteration_dirs(run_dir)
    f["iterations"] = its
    if iteration in ("final", None):
        chosen = its[-1] if its else None
    else:
        want = int(str(iteration).lstrip("0") or 0)
        chosen = next((p for n, p in its if n == want), None)
        if chosen is None and its:
            f["notes"].append(f"iteration {iteration} not present; using final")
            chosen = its[-1]
    f["iteration_dir"] = chosen[1] if isinstance(chosen, tuple) else chosen
    f["iteration_index"] = (its[-1][0] if iteration in ("final", None) and its
                            else (int(str(iteration).lstrip("0") or 0) if its else None))

    after = None
    if iteration not in ("final", None) and f["iteration_dir"]:
        after = _iteration_eqdsk(f["iteration_dir"])
    if after is None:
        hits = sorted(p for p in glob.glob(os.path.join(run_dir, "EQDSK*.OUT"))
                      if os.path.isfile(p))
        after = hits[0] if hits else None
    if after is None and f["iteration_dir"]:
        after = _iteration_eqdsk(f["iteration_dir"])
    f["eqdsk_after"] = after

    for name, cands in (
        ("baseline_dir", [cfg.get("baseline_dir"),
                          os.path.join(run_dir, "cheasebs_baseline"),
                          os.path.join(run_dir, "baseline"),
                          os.path.join(os.path.dirname(run_dir), "baseline")]),
        ("convergence_summary", [os.path.join(run_dir, "convergence_summary.json"),
                                 os.path.join(run_dir, "runs", "convergence_summary.json")]),
        ("acceptance", [os.path.join(run_dir, "cheasebs_acceptance.json")]),
        ("iteration_log", [os.path.join(run_dir, "iteration_log.csv"),
                           os.path.join(run_dir, "runs", "iteration_log.csv")]),
    ):
        f[name] = next((p for p in cands if p and os.path.exists(p)), None)

    f["summary"] = _json_or_empty(f["convergence_summary"])
    f["acceptance_rec"] = _json_or_empty(f["acceptance"])
    f["baseline_csv"] = None
    if f["baseline_dir"]:
        p = os.path.join(f["baseline_dir"], "profiles.csv")
        f["baseline_csv"] = p if os.path.isfile(p) else None
        p = os.path.join(f["baseline_dir"], "istar_current_profile_rhop.csv")
        f["baseline_istar_csv"] = p if os.path.isfile(p) else None
    else:
        f["baseline_istar_csv"] = None

    f["case"] = _case_from_dirname(os.path.basename(run_dir))
    return f


def _json_or_empty(path):
    if not path or not os.path.isfile(path):
        return {}
    try:
        return json.load(open(path))
    except Exception:
        return {}


def _case_from_dirname(name):
    """Axis and scale from the campaign's own directory naming, e.g.
    ``ne_ped_scale_1.300``. Used only to know which way a quantity should move;
    a directory that does not match simply gets no direction expectation."""
    m = re.match(r"^(?P<axis>.+?)_(?P<scale>\d+\.\d+)$", name)
    if not m:
        return {"axis": None, "scale": None}
    return {"axis": m.group("axis"), "scale": float(m.group("scale"))}


# --------------------------------------------------------------------------
# readers
# --------------------------------------------------------------------------

def read_gene_profile(path):
    """(rho_tor, T[keV], n[1e19 m^-3]) from a GENE profiles_<spec> file.

    Header wording differs between writers (TPED emits 'rho_tor', cheaseBS's
    bundled files 'rhot'), so only the column ORDER is relied on: radius,
    second radius, temperature, density.
    """
    data = np.loadtxt(path, comments="#")
    if data.ndim != 2 or data.shape[1] < 4:
        raise ValueError(f"{path}: expected 4 columns, got shape {data.shape}")
    return data[:, 0], data[:, 2], data[:, 3]


def thermal_pressure(profile_paths):
    """(rho_tor, p_th[Pa]) summed over whatever species are present.

    Species are interpolated onto the electron grid; a missing species is
    skipped rather than assumed zero-pressure silently -- the count comes back
    so the caller can say what went into the sum.
    """
    if not profile_paths.get("e"):
        return None, None, []
    rho, te, ne = read_gene_profile(profile_paths["e"])
    p = ne * te * PA_PER_1E19_KEV
    used = ["e"]
    for spec in ("i", "z"):
        path = profile_paths.get(spec)
        if not path:
            continue
        r_s, t_s, n_s = read_gene_profile(path)
        o = np.argsort(r_s)
        p = p + (np.interp(rho, r_s[o], n_s[o]) * np.interp(rho, r_s[o], t_s[o])
                 * PA_PER_1E19_KEV)
        used.append(spec)
    return rho, p, used


def read_csv_columns(path):
    """CSV with a header row into {column: float array}. Non-numeric -> NaN."""
    import csv as _csv
    with open(path, newline="") as fh:
        rows = list(_csv.reader(fh))
    if not rows:
        return {}
    head, body = rows[0], rows[1:]
    out = {k: np.empty(len(body)) for k in head}
    for i, row in enumerate(body):
        for k, v in zip(head, row):
            try:
                out[k][i] = float(v)
            except (TypeError, ValueError):
                out[k][i] = np.nan
    return out


def read_xy(path):
    """Two-column '# label label' text file (bootstrap diagnostic, I* target, q)."""
    d = np.loadtxt(path, comments="#")
    d = np.atleast_2d(d)
    return d[:, 0], d[:, 1]


def header_current(path):
    """Ip from EQDSK header line 4, field 1 -- the number the file itself claims."""
    with open(path) as fh:
        lines = [next(fh) for _ in range(5)]
    return float(lines[3][0:16])


# --------------------------------------------------------------------------
# derived equilibrium quantities
# --------------------------------------------------------------------------

class Equilibrium:
    """One EQDSK with the volume-weighted quantities the table needs.

    Geometry is done with the flux-surface integrals in TPED's
    equilibrium_geometry (Pappus volume off closed contours, agrees with
    CHEASE's ogyropsi to a fraction of a percent) rather than a cylindrical
    approximation, which is wrong on a spherical tokamak by a factor that grows
    with radius.
    """

    def __init__(self, path, want_geometry=True):
        self.path = path
        self.ds = _gfile_reader()(path).gfile_to_xarray()
        self.geom = None
        self.geom_error = None
        if want_geometry:
            try:
                self.geom = _geometry()(self.ds)
            except Exception as exc:                      # noqa: BLE001
                self.geom_error = f"{type(exc).__name__}: {exc}"

    # --- scalars off the header ---
    @property
    def r_axis(self):
        return float(self.ds.attrs["rmag"])

    @property
    def z_axis(self):
        return float(self.ds.attrs["zmag"])

    @property
    def b0(self):
        return abs(float(self.ds.attrs["Bctr"]))

    @property
    def ip(self):
        return header_current(self.path)

    # --- profiles ---
    @property
    def q(self):
        return np.asarray(self.ds["q"].values, dtype=float)

    @property
    def p(self):
        return np.asarray(self.ds["p"].values, dtype=float)

    @property
    def rho_tor(self):
        return np.asarray(self.ds.coords["rho_tor"].values, dtype=float)

    @property
    def rho_pol(self):
        return np.asarray(self.ds.coords["rho_pol"].values, dtype=float)

    def q_at(self, rho, coord="rho_tor"):
        x = self.rho_tor if coord == "rho_tor" else self.rho_pol
        o = np.argsort(x)
        return float(np.interp(rho, x[o], self.q[o]))

    def boundary(self):
        if "RBDRY" not in self.ds:
            return None
        r = np.asarray(self.ds["RBDRY"].values, dtype=float)
        z = np.asarray(self.ds["ZBDRY"].values, dtype=float)
        return (r, z) if r.size > 3 else None

    # --- volume-weighted ---
    def _vol_mask(self):
        v = np.asarray(self.geom["volume"].values, dtype=float)
        return v, np.isfinite(v)

    @property
    def volume(self):
        if self.geom is None:
            return np.nan
        v, m = self._vol_mask()
        return float(np.nanmax(v[m])) if m.any() else np.nan

    @property
    def a_minor(self):
        if self.geom is None:
            return np.nan
        return float(self.geom.attrs.get("a_minor", np.nan))

    @property
    def w_thermal(self):
        """3/2 int p dV over the whole plasma, from the EQDSK's own p(psi).

        This is total pressure as CHEASE stored it -- thermal plus the fixed
        fast-ion part -- so it is the stored-energy *response*, not a thermal
        decomposition. The thermal/fast split is checked separately, off the
        profiles and the baseline decomposition.
        """
        if self.geom is None:
            return np.nan
        v, m = self._vol_mask()
        if m.sum() < 5:
            return np.nan
        return float(1.5 * _trapz(self.p[m], v[m]))

    @property
    def p_volume_avg(self):
        if self.geom is None:
            return np.nan
        v, m = self._vol_mask()
        if m.sum() < 5 or not np.isfinite(self.volume) or self.volume <= 0:
            return np.nan
        return float(_trapz(self.p[m], v[m]) / (v[m].max() - v[m].min()))

    @property
    def beta_t(self):
        pav = self.p_volume_avg
        if not np.isfinite(pav):
            return np.nan
        return float(2.0 * MU0 * pav / self.b0 ** 2)

    @property
    def beta_n(self):
        """beta_t[%] * a[m] * B0[T] / Ip[MA]."""
        if not np.isfinite(self.beta_t) or not np.isfinite(self.a_minor):
            return np.nan
        ip_ma = abs(self.ip) / 1e6
        if ip_ma <= 0:
            return np.nan
        return float(self.beta_t * 100.0 * self.a_minor * self.b0 / ip_ma)

    def shape_at(self, psi_norm=0.95):
        """(kappa, delta, R_geo, r_minor) on the psi_N = 0.95 surface."""
        if self.geom is None:
            return {}
        pn = np.asarray(self.geom.coords["psi_norm"].values, dtype=float)
        out = {}
        for key, name in (("kappa", "kappa"), ("delta_upper", "delta_u"),
                          ("delta_lower", "delta_l"), ("R_geo", "R_geo"),
                          ("r_minor", "r_minor")):
            v = np.asarray(self.geom[key].values, dtype=float)
            m = np.isfinite(v) & np.isfinite(pn)
            out[name] = (float(np.interp(psi_norm, pn[m], v[m])) if m.sum() > 3
                         else np.nan)
        if np.isfinite(out.get("delta_u", np.nan)):
            out["delta"] = 0.5 * (out["delta_u"] + out["delta_l"])
        return out

    def r_geo_at(self, rho_tor=0.5):
        """Flux-surface geometric centre at one radius -- the Shafranov-shift probe."""
        if self.geom is None:
            return np.nan
        x = np.asarray(self.geom.coords["rho_tor"].values, dtype=float)
        v = np.asarray(self.geom["R_geo"].values, dtype=float)
        m = np.isfinite(v) & np.isfinite(x)
        if m.sum() < 4:
            return np.nan
        o = np.argsort(x[m])
        return float(np.interp(rho_tor, x[m][o], v[m][o]))


def boundary_displacement(eq_a, eq_b, n_theta=512):
    """RMS and max |dr| between two LCFS contours, sampled at common angle.

    Both boundaries are re-parameterised by poloidal angle about the *source*
    magnetic axis, because the two files do not share a point count (the source
    EFIT and CHEASE's 513x513 output rarely do) and comparing index-by-index
    would measure the resampling, not the boundary.
    """
    ba, bb = eq_a.boundary(), eq_b.boundary()
    if ba is None or bb is None:
        return None
    r0, z0 = eq_a.r_axis, eq_a.z_axis
    grid = np.linspace(-np.pi, np.pi, n_theta, endpoint=False)

    def sample(rz):
        r, z = rz
        th = np.arctan2(z - z0, r - r0)
        rad = np.hypot(r - r0, z - z0)
        o = np.argsort(th)
        th, rad = th[o], rad[o]
        # wrap one period either side so the interpolation is periodic
        th_ext = np.concatenate([th - 2 * np.pi, th, th + 2 * np.pi])
        rad_ext = np.concatenate([rad, rad, rad])
        return np.interp(grid, th_ext, rad_ext)

    ra, rb = sample(ba), sample(bb)
    d = np.abs(ra - rb)
    return {"rms_m": float(np.sqrt(np.mean(d ** 2))), "max_m": float(d.max())}


# --------------------------------------------------------------------------
# check plumbing
# --------------------------------------------------------------------------

def row(name, group, expect, reference, value, delta=None, tol=None,
        verdict=SKIP, note="", fmt="{:.4g}"):
    return {"check": name, "group": group, "expect": expect,
            "reference": reference, "value": value, "delta": delta,
            "tol": tol, "verdict": verdict, "note": note, "fmt": fmt}


def _rel(new, ref):
    if ref in (None, 0) or new is None:
        return np.nan
    if not (np.isfinite(new) and np.isfinite(ref)) or ref == 0:
        return np.nan
    return (new - ref) / abs(ref)


def _grade(value, warn, fail):
    """Absolute-value grading for an invariant: <= warn PASS, <= fail FLAG."""
    if value is None or not np.isfinite(value):
        return SKIP
    v = abs(value)
    if v <= warn:
        return PASS
    return FLAG if v <= fail else FAIL


def _direction(delta, scale, dead=1e-3, gated=True):
    """Did a quantity move, and the way the scan asked it to?

    dead is the "it never moved" floor: on a run whose profiles really were
    scaled, a relative change under 0.1% is the signature of profiles that
    never reached CHEASE, which is a FAIL rather than a small effect.

    gated=False is for identity replays and for runs where no profile change
    could be measured -- there, nothing is *supposed* to move, so the row is
    reported instead of scored. Passing the scan's own measured p_th change in
    as the gate is what keeps "unchanged" from meaning two opposite things.
    """
    if delta is None or not np.isfinite(delta):
        return SKIP, "no value"
    if not gated:
        return PASS, "reported, not gated: no profile change measured for this run"
    if abs(delta) < dead:
        return FAIL, (f"unchanged (|d| < {dead:.3g}) -- the scaled profiles may "
                      "never have reached CHEASE")
    if scale is None:
        return PASS, "moved; no scan direction known for this directory"
    want_up = scale > 1.0
    got_up = delta > 0
    if want_up == got_up:
        return PASS, ""
    return FAIL, f"moved the wrong way for scale {scale:g}"
