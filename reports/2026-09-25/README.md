# NSTX CHEASE-BS scaling review - 2026-09-25

| Discharge | PowerPoint | PDF (all slides) |
|---|---|---|
| 129015 | [26 slides](NSTX-129015-2026-09-25.pptx) | [26 pages](NSTX-129015-2026-09-25.pdf) |
| 129038 | [16 slides](NSTX-129038-2026-09-25.pptx) | [16 pages](NSTX-129038-2026-09-25.pdf) |

PDFs preserve the reviewed slide renderings as full-page images, with slide bookmarks. PowerPoints retain the original presentation. Both separate mtanh and apply_omX, use two slides per available variant, and finish each family with a pressure-chain comparison. Colors identify the same scaling type/factor throughout; comparison borders identify solver settings.

## Scope

Saved data: `../../data/129015_129038/pm0.2_ncscal{1,4}_warmup{0,1,2}`. Grid: 0.80, 0.90, 0.95, 1.00, 1.05, 1.10, 1.20, with Te/ne mtanh and omt/omne transforms. There are 225 saved convergence summaries out of 336 intended cases. Missing summaries are unavailable, not counted as failed runs. Particularly sparse: NCSCAL=1 on 129038 (one saved result).

Saved settings match on I-star replay, bootstrap/I-star mixing 0.1/0.05, Ip tolerance 0.2%, bootstrap/q tolerances 1%, amplitude tolerance 0.0001, maximum 50 iterations, and regularization. NCSCAL is verified from campaign namelists; warmup from saved configs. Matching settings do not establish executable build identity.

## 129015

- **NCSCAL=4:** all 84 cases (28 per warmup) report convergence in exactly two iterations, 00/01. Warmup=0 makes only small amplitude changes; warmup=1/2 retain amplitude 1. This is repeatable early stopping, not an independent demonstration of closure.
- **NCSCAL=1:** warmup=0 gives 17/22 reported convergences, only one in two iterations; final amplitudes span 0.030-2.532. Warmup=1/2 each give 3/18 reported convergences, with amplitude frozen at 1. Several bootstrap/q traces decrease while Ip remains above tolerance. With warmup=0, some errors initially improve then grow or oscillate; more iterations do not always mean progress.
- **Pressure chains:** NCSCAL=4 produces compact, largely overlapping warmup curves, but the source baseline lies above the reconstructed q chains and below the shear chains. The offset remains near scale=1. In mtanh, Te and ne can give opposite q slopes; steep-point shear has turns. In apply_omX, q generally decreases with pressure, while shear decreases at T_top and increases at n_top; the steep point has more curvature. No universal q/shear slope applies.
- **Outliers:** failed NCSCAL=1/warmup=0 outputs reach very large solved/source pressure ratios and irregular q/shear excursions, stretching the combined plots. Their crosses mark nonconvergence; the connecting lines do not make them trustworthy scan points.

## 129038

- **NCSCAL=4:** warmup=0 gives 28/28 reported convergences, 25 in two iterations. Warmup=1/2 each give 25/27, 24 in two iterations. Their failed ne x0.8 and omt x1.2 runs reach 50 iterations; Te x0.8 is unavailable for those two warmups.
- **Error trends:** the failed ne x0.8 Ip error rises then declines to about 0.274%, still above the 0.2% tolerance. The failed omt x1.2 trace is oscillatory and ends at about 0.895%. A late downward trend therefore does not imply convergence. Warmup=0 permits amplitude adjustment and reports convergence for these cases.
- **Pressure chains:** most NCSCAL=4 warmup curves overlap closely; difficult edge cases separate. The source baseline is close to the n_top chains but offset from T_top/steep shear chains. For mtanh, q generally falls with pressure at n_top but rises at T_top/steep. For apply_omX, q generally falls at all three points; shear falls at n_top and rises at T_top/steep. The nonconverged omt x1.2 point turns upward in q, disrupting the otherwise decreasing trend.
- **NCSCAL=1:** only Te x0.8/warmup=0 is available (23 iterations, final amplitude 0.0478). It cannot establish a normalization-wide trend.

## Interpretation

The source star is the original equilibrium. The hollow diamond is fitted baseline pressure at source q/shear, not a reconstructed scale=1 equilibrium. Pressure chains show solved/source pressure at fixed reference Chapman radii, not solved/requested pressure error. Compare each family with its own reconstructed unity case, and evaluate pressure adherence and closure separately before accepting profiles for physics scans.

The decks use TPED GFileData and shear helpers. Chapman radii: 129015 T_top=0.8425, n_top=0.8745, steep=0.9515; 129038 n_top=0.6737, T_top=0.8447, steep=0.8535. The pressure/q source arrays agree within each discharge. No new solver runs or solver code changes were made for this report.
