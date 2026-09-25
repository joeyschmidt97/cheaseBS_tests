# NSTX 129015 CHEASE-BS review — 2026-09-24

- [PowerPoint](NSTX-129015-review.pptx): final 14-slide review.
- [Slide preview](index.html): open locally in a browser; no server required.
- [Brief advisor update](advisor-update.md).
- [Run provenance and sampled pressure-chain data](data/run_manifest.json).
- [Fixed scaling colors](data/scaling_colors.json).

## Layout and coverage

Slides 1–9: mtanh pedestal-height scans (Te and ne), four NCSCAL/warmup combinations. Slides 10–14: apply_omX (omt and omne), two NCSCAL=4 combinations. Each variant gets transformations/error traces followed by iteration bars/pressure chains. Each family ends with an overlay using marker borders for the solver variant. No NCSCAL=1 apply_omX data exist in the selected campaigns.

46 saved runs: runs_smallscale/129015_20260901_13-45-25 (8), runs_amptest/129015_20260903_15-51-42 (7), data/pm0.3_ncscal4/129015 (12), data/pm0.3_ncscal4_warmup0/129015 (19). The first two campaigns are in the sibling NSTXU_lithium_study repository under reshape_convergence. Nested copies and separate _omit attempts are excluded. No new equilibrium solves or solver modifications were performed.

## Plot conventions

One fixed color for each scaling type/factor is used in transformation curves, iteration traces, bars and pressure-chain marker fills, across solver variants. Marker shapes identify scaling type. Marker borders and connecting lines on the final overlays identify solver variant. Black profile curves are the reference; solid profiles are electrons, dashed profiles ions where shown. The plots use the saved profiles_e/profiles_i data directly. Older campaign panels show electrons, as in their original comparison images; NCSCAL=4 panels retain electron and ion profiles.

Iteration counts include 00. Hatched bars and crosses identify solver-reported nonconvergence. Missing first-iteration bootstrap/q changes are not fabricated. Error curves are relative fractions, with a small linear region below 1e-7 to preserve zeros. All-zero traces use linear axes. Amplitude change may be explicitly set to zero when the solver declares convergence; it is not evidence of a satisfied independent amplitude gate.

Pressure chains use solved total pressure divided by source pressure at common Chapman rho_tor values: T_top=0.8425, n_top=0.8745, steep=0.9515. Points were refitted once from the shared reference electron profile using TPED chapman_points.py. EQDSKs and shear use GFileData and cheasebs_runner.py helpers. The star is the source equilibrium; the hollow diamond is the fitted baseline pressure evaluated at source q/shear, not a separate reconstructed equilibrium. Curves join pressure-ordered points within one scaling family and solver variant. No cross-family slope is fitted. All source equilibria agree in pressure and q.

## Limits

NCSCAL=1 follows the archived config's reference to chease_namelist_nstx and that preset's current contents; its original runtime namelist is unavailable. NCSCAL=4 is explicit in the saved campaign namelist. Older campaigns use I-star mixing 0.02, Ip tolerance 2%, and bootstrap/q tolerance 0.1%; newer campaigns use 0.05, 0.2%, and 1%, respectively. Scale coverage also differs. These plots do not isolate the causal effect of NCSCAL or certify equilibrium validity.

Manifest source-path roots ${CHEASEBS_TESTS} and ${NSTXU_LITHIUM_STUDY} denote the corresponding repository roots; remaining path components preserve their original Windows notation. Slide notes retain the original source paths. This folder is a portable review snapshot, not a solver or plotting environment.

The broader 129038 example in the advisor update is outside this deck: data/pm0.3_ncscal4_warmup0/129038/omt_1.300. It reaches 50 iterations without convergence. Pressure relative L2 error is 42.2019%, solved/requested pressure span is 0.586486, and solved/source flux span is 0.586616. Target pressure = active thermal + max(source pressure − reference thermal, 1 Pa), evaluated on the source rho_pol grid; the L2 is sample-weighted, not volume-weighted.

## Validation

All 46 manifest run records and pressure-chain samples match the preceding review. 936 plotted artist colors were checked against the 24-entry fixed palette. The 14-slide PPTX passed package and layout validation; rendered slides were visually reviewed. This does not constitute execution or testing in Microsoft PowerPoint.
