# Advisor update: NSTX CHEASE-BS scaling tests

I compared the saved NSTX 129015 runs using I-star replay:

- **NCSCAL=1, warmup=2:** 7/8 runs stop at 00/01. Amplitude stays at 1; ne pedestal ×0.9 reaches 50 iterations without convergence.
- **NCSCAL=1, warmup=0:** amplitude can change. All 7 report convergence (5 in two iterations); Te pedestal ×0.7 takes 6 iterations and finishes at amplitude ≈0.20.
- **NCSCAL=4, warmup=2:** 11/12 stop at 00/01 with amplitude fixed at 1; ne pedestal ×0.7 stalls at 50 iterations.
- **NCSCAL=4, warmup=0:** all 19 stop at 00/01 despite allowing amplitude changes (final amplitudes ≈0.9996–1.0009). Removing warmup does not remove the early exits.

Other observations:

- With NCSCAL=4, the source baseline sits off the reconstructed pressure–q trend; shear is also offset. Some chains turn rather than remaining monotonic. Source baseline and reconstructed scale=1 must be compared separately.
- Early stopping does not establish pressure adherence or full self-consistency: the current stopping rule checks Ip and successive bootstrap/q changes, without a pressure-error gate. Across the broader NSTX tests, one failed 129038 omt ×1.3 run has ≈42% pressure-profile error, with pressure-span loss tracking flux-span loss.
- These are campaign comparisons, not isolated NCSCAL tests: mixing, tolerances and sampled scales differ. NCSCAL=1 is inferred from the referenced preset. Next is a matched-settings test retaining I-star, with pressure-adherence and convergence checks.
