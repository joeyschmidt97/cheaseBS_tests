# cheaseBS scale-mapping and response-symmetry test

One symmetric sweep answering two questions left open by the 2026-08-18
five-point 132543 sparse-grid run. Both are properties of the same curve, so
they share a sweep rather than needing two.

## Q1 — what does a scale factor actually do?

`apply_mtanh_ped(scale_height=s)` multiplies the mtanh **step amplitude**
`a_y0`, not the pedestal-top value, which is `y_sep + a_y0 + core`. Measured on
the 132543 seed:

| axis setting | measured pedestal-top change |
|---|---|
| `ne_ped_scale = 1.3` | **+14.4%** |
| `ne_ped_scale = 0.7` | **-14.4%** |
| `Te_ped_scale = 1.3` | **+21.3%** |

So a scan box quoted in scale factors is far narrower than the same numbers
read as pedestal-top fractions, and the survey's absolute bounds table
(`T_e,ped` 0.2-0.8 keV, `n_e,ped` 3-7e19 m^-3) cannot be converted into scale
factors without this curve. The sweep produces a slope per variable; a constant
slope column means the mapping is linear and one number per variable does the
conversion.

## Q2 — why is the response asymmetric?

| point | iterations | Ip error |
|---|---|---|
| `ne x1.3` | **19** | 1.62% -> **0.77%** |
| `ne x0.7` | **2** | 1.62% -> 1.62% (unmoved) |

The profile input is exactly symmetric (+/-14.4% at the pedestal top), so the
asymmetry is in cheaseBS's response. Candidates:

- the floored `p_fast` split (`p_total = p_thermal(new) + p_fast(reference, floored)`)
- quasineutrality clipping `nz = (ne - ni)/qz` at low density
- `qspec=on` pinning the current asymmetrically

## Running it

```bash
python run_scale_symmetry_test.py --dirpath /path/to/132543
```

Options: `--scales 0.7,0.85,1.0,1.15,1.3` (must stay symmetric about 1.0),
`--vars ne` to sweep one variable, `--outroot` for the output location.

Default is 5 scales x 2 variables = 10 cheaseBS runs at up to 25 iterations
each. Reduce with `--vars ne` first if cluster time is tight: Q2 is a density
question and Q1 needs only one variable to establish whether the mapping is
linear at all.

## Reading the output

`runs/<timestamp>/scale_symmetry_results.json` holds every row. Per point,
`runs/<timestamp>/<tag>/` has the equilibrium, the acceptance record,
`iteration_log.csv` and `iteration_errors.png`.

The two printed tables are the answers. For Q1, look at whether the `slope`
column is constant. For Q2, compare each `+delta` row against its `-delta`
partner: equal iteration counts and equal-and-opposite Ip shifts mean the
response is symmetric and the 19-vs-2 split was specific to that one point; a
systematic difference points at a one-sided mechanism.

## Requirements

NERSC, with the compiled CHEASE binary reachable through
`TPED/config/user_config.yaml`, and TPED on `gene_datatree` at or past
`036087d` (the run needs `plot_errors` and the fixed iteration-log column
matching).
