# Superseded — moved to `NSTXU_lithium_study/reshape_convergence/`

**Do not run `run_reshape_convergence.py` in this directory.** It drifted away from
the campaign and would measure a reconstruction nobody will run:

- it calls `apply_mtanh_ped` (pedestal-only) while the campaign froze on
  `apply_mtanh_full` (Stefanikova full-profile) on 2026-08-23;
- it hardcodes the discharge list, data root, pedestal window, fit settings and
  scan scales, all of which `NSTXU_lithium_study/pedestal_scan.py` now owns as
  the single source of truth;
- its 129038 analysis radius is flagged `placeholder_radii` and its verdicts on
  that shot were never meaningful.

The replacement imports `pedestal_scan` instead of restating it, so the same
drift cannot happen again.

## `runs/20260820_10-34-35/` contains no convergence data

Every row in `reshape_convergence.json` carries

    "error": "ValueError: Incorrect config path key used \"CHEASE_PATH\" ..."

with `wall_s` of 1–40 ms. cheaseBS never ran; this is a harness check that was
left on disk looking like a result set. The directory layout, the EQDSK-named
subdirectories and the populated JSON all suggest otherwise — read the `error`
key before drawing anything from it.

The replacement prints an explicit `NO CONVERGENCE DATA` banner when every
solve raises, so this failure mode cannot be mistaken for a result again.
