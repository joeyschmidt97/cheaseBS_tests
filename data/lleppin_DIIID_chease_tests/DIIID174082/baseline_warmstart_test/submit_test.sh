#!/bin/bash
#SBATCH --job-name=warmstart_test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --constraint=cpu
#SBATCH --qos=debug
#SBATCH --account=m2116
#SBATCH --output=log.warmstart_test.%j.out
#SBATCH --error=log.warmstart_test.%j.err

set -e
module load python
source ~/.bashrc
conda activate lleppinenv

cd /global/homes/l/lleppin/TPED/projects/cheaseBS
CFGDIR=/global/homes/l/lleppin/DIIID174082/baseline_warmstart_test

for c in omt0p7_omne0p7 omt0p8_omne1p0 omt1p1_omne1p1; do
  echo "=== $c ==="
  python run_chease_iterative_profiles.py --config "$CFGDIR/baseline_warmstart_$c.json" \
    > "$CFGDIR/baseline_warmstart_$c.log" 2>&1 || echo "$c: FAILED (exit $?)"
  echo "=== $c done ==="
done
