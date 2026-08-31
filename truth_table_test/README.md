# cheaseBS profile-scaling truth table

A reader, not a runner. Point it at cheaseBS runs that already exist and it says,
row by row, whether the quantities that were **imposed** survived the profile
scan and whether the quantities the scan was supposed to **drive** actually
moved -- and moved the right way.

The two groups fail in opposite directions, which is the whole reason this
exists. A control that changed means the comparison is not a profile scan. A
response that did *not* change means the scaled profiles never reached CHEASE --
the 2026-08-26 reference-profile bug, which read as a converged run for weeks.
One "converged" flag cannot separate those.

## Running it

```bash
# one point
python run_truth_table.py --run /path/to/Te_ped_scale_1.300

# one campaign, plus its cross-case ordering checks
python run_truth_table.py --campaign /path/to/reshape_convergence/runs/129015_20260824_20-28-57

# a whole tree: one directory per shot per timestamp, four box-edge points each
python run_truth_table.py --campaign /path/to/reshape_convergence/runs --summary-only

# score the documented NSTX replay point rather than the final iteration
python run_truth_table.py --campaign <dir> --iteration 00
```

`--campaign` is **recursive**. Any directory carrying an `EQDSK*.OUT`,
`iteration_NN/` or `cheasebs_run_config.json` is one run point; the walk does not
descend into a point's own artifact subdirectories. Points are grouped by parent
directory, so a tree of four shots gets one ordering block per campaign rather
than one meaningless block spanning all of them. `--summary-only` collapses each
point to a line and keeps the cross-case blocks -- the full tables run to pages
at sixteen points.

Options: `--source-eqdsk` when the source copy is missing from the run dir,
`--tol-ip-rel` and `--dead` to override the two tolerances worth overriding,
`--quiet-notes` for the table without explanations, `--out` for the JSON.

Exit status is 0 when nothing FAILed, 1 otherwise, so a wrapper can tell "ran"
from "verified". A JSON copy of every row goes to `truth_table.json` in the
campaign root (or beside a single `--run`).

`TPED` must be importable: it supplies the EQDSK reader and the flux-surface
geometry. Set `TPED_ROOT` to the directory *containing* the checkout if a bare
`import TPED` does not work.

### Which iteration

`--iteration` chooses only which equilibrium the **control** and **response**
blocks are scored on; `final` is the default. The convergence block always reads
the whole path regardless.

Iteration 00 is worth scoring separately because cheaseBS's own documentation
names it the validated NSTX `I*` replay point and says later outer iterations
are numerically stable but drift away from it. So run both and quote which one
you mean; `iter00_vs_final_drift` puts a number on the gap instead of leaving it
as a caveat.

## What it reads

Nothing has to be re-derived, and no separate baseline solve is needed -- a
cheaseBS run directory already carries both sides of its own comparison:

| Artifact | Used for |
|---|---|
| `g<shot>.<time>` | the frozen **source** equilibrium (the "before") |
| `EQDSK*.OUT`, or `iteration_NN/artifacts/EQDSK*.OUT` | the reconstruction (the "after") |
| `profiles_{e,i,z}` | active thermal pressure |
| `reference_profiles/REF_profiles_{e,i,z}` | reference thermal pressure |
| `cheasebs_baseline/profiles.csv` | reference `p_fast`, `p_th`, bootstrap current |
| `cheasebs_baseline/istar_current_profile_rhop.csv` | reference total / component `I*` |
| `iteration_NN/bootstrap_diagnostic_rhop.txt` | this iteration's bootstrap current |
| `iteration_NN/istar_target_rhop.txt` | the total `I*` actually replayed |
| `convergence_summary.json`, `cheasebs_acceptance.json` | targets, iterations, q errors |

`--campaign` also reads `reshape_convergence.json` when present, for the axis
and scale of each point -- that is what tells the checker which way each
response is *supposed* to move. Points whose recorded `savedir` is on another
machine are reported as missing rather than silently skipped.

## The table

Verdicts: **PASS** met, **FLAG** needs a human, **FAIL** breaks the scan,
**SKIP** artifact absent. FLAG never sets the exit status.

### Controls -- must not move

| Row | Reference | Criterion |
|---|---|---|
| `run_completed` | -- | an `EQDSK*.OUT` exists; FAIL if the point produced no equilibrium |
| `solver_state` | `max_iter` | FLAG when the run hit the iteration cap, or when cheaseBS's own `converged` is False (the acceptance policy does not gate on that flag, but "slow" and "diverging" must be distinguishable) |
| `ip_closure` | `target_ip_a` | \|dIp\|/Ip <= the run's own `tol_ip_rel` (else 0.02) |
| `bt_invariant` | source \|Bctr\| | identical to 1e-6; a profile scan must not touch the imposed field |
| `lcfs_rms/a` | source LCFS | rms displacement <= 0.005 a (FLAG to 0.02). **The fixed-boundary gate.** Both boundaries are resampled at common poloidal angle about the source axis, because the source EFIT and CHEASE's 513x513 output do not share a point count |
| `lcfs_max/a` | source LCFS | report/FLAG only -- the X-point region dominates the maximum, and a hard bound there rejects sound reconstructions |
| `kappa_95` | source | relative change <= 0.02 |
| `delta_95` | source | absolute change <= 0.02, mean of upper and lower |
| `plasma_volume` | source | relative change <= 0.02. That bound is set by the grid change alone: two identity replays (132588, 129015) move the volume 0.8-1.2% with the boundary held fixed |
| `pfast_reference_frozen` | `REF_profiles` | FAIL when the reference profiles are aliased onto the run profiles -- then the baseline was rebuilt from the scaled profiles, `p_fast` moved with the scan, and this is not a thermal-only scan |
| `istar_baseline_closure` | baseline CSV | `max|reconstructed - total| / max|total|` <= 1e-9 (documented as machine precision) |
| `physicality` | -- | q finite and nonzero, p >= 0, flux surfaces contourable |

### Responses -- must move, and the right way

Scored against the scan direction taken from the run directory name or the
campaign JSON. On an identity replay, or when no profile change can be measured,
these rows report instead of scoring -- "unchanged" must not mean two opposite
things.

| Row | Reference | Criterion |
|---|---|---|
| `p_th_profile_scale` | `REF_profiles` (or baseline `p_th_pa`) | median `p_th` ratio over rho_tor <= 0.995. Must move by > 0.1% and in the scan's direction. **The ingestion test with teeth** |
| `p_axis_ingestion` | `p_fast(reference) + p_th(active)` | axis pressure of the reconstruction within 5% of that sum (FLAG to 15%). Direct check that CHEASE solved with the profiles you think it did |
| `p_axis_response` | source p(0) | moved, in the scan's direction |
| `beta_t` | source | `2 mu0 <p>_V / B0^2`, volume-averaged over the real flux-surface volumes; moved in the scan's direction |
| `beta_N` | source | report only |
| `W_stored` | source | `3/2 int p dV` on the EQDSK's own total p -- thermal plus the fixed fast-ion part; moved in the scan's direction |
| `axis_R_shift_mm` | source `rmaxis` | must exceed a 0.1 mm grid-noise floor to count as a resolved response |
| `axis_Z_shift_mm` | source `zmaxis` | report only; large vertical motion at fixed boundary is suspect |
| `shafranov_rho0.5_mm` | source | flux-surface centre at rho_tor 0.5, outward at higher pressure |
| `q_at_analysis_radii` | source q | recorded, gated only as a blowup bound (0.5). A reshaped pedestal is *supposed* to move q |
| `q_edge_error` | source / QSPEC | 1% is the documented figure; these NSTX cases carry a standing 10-20% offset on the outermost surface (grid/X-point artifact), so >1% is a FLAG to read |
| `q0_qmin_q95` | source | report only; read for discontinuities or a rational-surface crossing |
| `bootstrap_response` | baseline decomposition | `int |j_bs| drhop` moved, in the scan's direction. A shape proxy, not `I_bs` in amps |
| `istar_total_response` | baseline total `I*` | `rms(I*_target - I*_baseline)/max|I*_baseline|` > 1e-3. Zero means the replayed current never saw the profile change |

### Convergence -- how the loop got there

Read from `iteration_log.csv` and the per-iteration EQDSKs, against the run's own
configured tolerances rather than this checker's guesses. Produced for every
`--iteration` choice.

| Row | Reference | Criterion |
|---|---|---|
| `iterations_run` | `max_iter` | FLAG at the cap: a residual still falling was cut off, and "slow" cannot be told from "diverging" |
| `ip_error` | `tol_ip_rel` | final `ip_error_rel` <= tol (FAIL past 5x) |
| `bootstrap_change` | `tol_bs` | final `bootstrap_change_rel` <= tol |
| `q_change` | `tol_q` | final `q_change_rel` <= tol |
| `amplitude_change` | `tol_a` | final `amplitude_change_rel` <= tol |
| (all four) | the trace | each row also reports first -> final, the saturation iteration, and ringing: >= 3 sign flips in the tail downgrades a passing row to FLAG, and a residual that grew from its first value is a FAIL. Saturation is not convergence -- a loop can saturate at a value that never meets its tolerance, which is what cheaseBS does with `tol_ip_rel` under `qspec=on` |
| `driver_converged_flag` | -- | report only. False with every residual passing means one unmet tolerance; False alongside a FLAGged residual is the real thing |
| `iter00_vs_final_drift` | iteration 00 | max relative difference in `Ip`, axis pressure and `q(0.985)` between the first and last iteration. <= 1% PASS, FLAG to 10%. This is the documented NSTX drift, measured -- large is a reason to quote one iteration explicitly, not a rejection |

### Across cases -- only a campaign can answer these

| Row | Criterion |
|---|---|
| `Bt target identical` | the imposed field is the same number in every case |
| `Ip target identical` | so is the target current; a target that moves between cases invalidates the comparison however well each point scores alone |
| `<axis>/<quantity>` | low scale below high scale in `p_th`, axis pressure, `beta_t`, `W_stored`, bootstrap. A single point can look plausible while the pair is ordered backwards |

## Reading a result

The first two questions, in order:

1. **Any FAIL in the controls?** Then stop -- the boundary, the field, the
   current target or the fast-pressure reference moved, and no physics should be
   read off the scan.
2. **Any FAIL in the responses?** A response that did not move is the profiles
   not reaching CHEASE. A response that moved the wrong way is the scale not
   being applied as intended, or another setting changing between points.

Only then are the FLAGs worth time: the edge-q offset, an iteration cap, a
response under the grid-noise floor.

Two specifics on these NSTX cases:

- **Iteration 00 is the documented replay benchmark** on the `I*` branch; later
  outer iterations are numerically stable but drift. Score both
  (`--iteration 00` and the default) and quote them separately.
- **A constrained q is not an unchanged equilibrium.** With QSPEC on and the
  current targeted, a working pressure scan can hold `Ip`, `Bt`, the LCFS and
  edge q nearly fixed while `p`, `beta`, the bootstrap current, the axis
  position and the internal flux geometry all move for real. That pattern is the
  expected outcome, not a null result.

## What it does not do

- No mesh-convergence check. "The effect exceeds the production-to-fine mesh
  difference" needs a second solve at higher resolution, which is a run, not a
  read. The 0.1 mm axis floor and the 0.02 volume bound are stand-ins measured
  from identity replays.
- `bootstrap_response` is a profile-shape integral, not `I_bs` in amps, so its
  magnitude is not a bootstrap fraction.
- `W_stored` is built from the EQDSK's total pressure. The thermal/fast split is
  checked separately, off the profiles and the baseline decomposition.
