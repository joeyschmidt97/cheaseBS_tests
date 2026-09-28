#!/bin/bash
#SBATCH --job-name=cheasebs_scan
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --time=01:00:00
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --account=m2116
#SBATCH --output=log.cheasebs_scan.%j.out
#SBATCH --error=log.cheasebs_scan.%j.err

# Usage: sbatch submit_scan.sh /path/to/scan_config.json
#
# --time/--account/--qos above are placeholders copied from astra_runner's
# run_astra.sh -- adjust to your allocation and to how many cases * ~10min
# the scan actually needs given n_workers in scan_config.json.
#
# --nodes=1 already gets the whole node exclusively under NERSC's CPU-node
# policy (no sharing on --constraint=cpu); --ntasks=1/--cpus-per-task=128
# don't change that, they just document that this is one Python process
# (run_scan.py) that internally spawns up to n_workers CHEASE subprocesses
# via a process pool, rather than something SLURM should split into tasks.
# n_workers in scan_config.json is the actual concurrency control -- CHEASE's
# peak RSS measured at ~1GB/process, so up to 128 concurrent runs fits
# comfortably in the node's 256GB even with margin for oversized cases.

set -e

module load python
source ~/.bashrc
conda activate lleppinenv

# Hardcoded rather than derived from $0/${BASH_SOURCE[0]}: SLURM copies the
# submitted script into a per-job spool directory before running it, so
# that would resolve to the spool dir instead of this script's real
# location on the shared filesystem.
RUN_SCAN=/global/homes/l/lleppin/TPED/projects/profile_tools/scripts/run_scan.py
python "$RUN_SCAN" --config "$1"
