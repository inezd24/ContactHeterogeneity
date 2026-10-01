#!/bin/bash

##########################################################################################################################
# Script by Inez Derkx, contact: inez.derkx@swisstph.ch
#
# run_pipeline.sh: submits the full network-generator analysis to SLURM as a chain of dependent jobs.
#
# Run this from a LOGIN node with plain bash (it only submits jobs, it does no computing itself):
#
#   bash run_pipeline.sh          # everything
#   bash run_pipeline.sh main     # 01 -> 01a -> 02 -> 03a / 03b / 04a / 04b
#   bash run_pipeline.sh sims     # 03a / 03b / 04a / 04b only (graphs from 02 already exist)
#   bash run_pipeline.sh other    # 05 -> 06 (additional locations)
#
# Dependency structure:
#
#   01 ──> 01a ──> 02 ──┬──> 03a  SEIR randomized
#                       ├──> 03b  SEIR factorial
#                       ├──> 04a  SIS randomized
#                       └──> 04b  SIS factorial
#
#   05 ──> 06            (runs in parallel with the main chain)
#
# NOTES:
# 1. 00_constructNetworksFunction.R is not submitted: it is sourced by the other R scripts.
# 2. 02 and 06 have no .sh file, so they are submitted here with --wrap. Adjust their memory/time below.
# 3. For array jobs, a dependency waits for ALL array tasks. With DEP_TYPE="afterok", one failed task
#    cancels everything downstream. Use "afterany" if downstream steps can cope with missing tasks.
# 4. Flags given here (--chdir, --dependency) override nothing in the .sh files except where named;
#    all #SBATCH settings (array size, memory, time) are still taken from each .sh file.
##########################################################################################################################

set -euo pipefail

### SETTINGS ###

PROJECT_DIR="/scicore/home/chitnis/derkx0000/Derkx2026_publication"
LOG_DIR="${PROJECT_DIR}/SLURM"
R_MODULE="R/4.2.1-foss-2022a"
DEP_TYPE="afterok"        # or "afterany"

# Resources for the two collection scripts that have no .sh file
COLLECT_MEM="16G"
COLLECT_TIME="02:00:00"

MODE="${1:-all}"

case "$MODE" in
  all|main|sims|other) ;;
  *) echo "Unknown mode '$MODE'. Use: all | main | sims | other"; exit 1 ;;
esac

cd "$PROJECT_DIR"
mkdir -p "$LOG_DIR"


### PRE-FLIGHT CHECKS ###

REQUIRED=(00_constructNetworksFunction.R
          01_NetworkGeneratorComparison.sh 01a_Run_NCRG.sh 02_collectGraphs.R
          03a_SEIR_Randomized.sh 03b_SEIR_Factorial.sh 04a_SIS_Randomized.sh 04b_SIS_Factorial.sh
          05_Graphs_Other_Locations.sh 06_collectGraphsOtherLocations.R)

for f in "${REQUIRED[@]}"; do
  [[ -f "$f" ]] || { echo "ERROR: missing $f in $PROJECT_DIR"; exit 1; }
done


### HELPERS ###

# Build the dependency flags from a colon-separated list of job IDs (empty = no dependency)
dep_flags() {
  local deps="${1:-}"
  if [[ -n "$deps" ]]; then
    echo "--dependency=${DEP_TYPE}:${deps} --kill-on-invalid-dep=yes"
  fi
}

# Submit an existing .sh file; prints the job ID
submit_sh() {
  local script="$1" deps="${2:-}"
  # shellcheck disable=SC2046
  sbatch --parsable --chdir="$PROJECT_DIR" $(dep_flags "$deps") "$script" | cut -d';' -f1
}

# Submit an R script that has no .sh file; prints the job ID
submit_r() {
  local rscript="$1" name="$2" deps="${3:-}"
  # shellcheck disable=SC2046
  sbatch --parsable \
         --chdir="$PROJECT_DIR" \
         --job-name="$name" \
         --output="${LOG_DIR}/%x.%j.out" \
         --error="${LOG_DIR}/%x.%j.err" \
         --cpus-per-task=1 \
         --mem="$COLLECT_MEM" \
         --time="$COLLECT_TIME" \
         $(dep_flags "$deps") \
         --wrap="ml purge; ml ${R_MODULE}; Rscript ${rscript}" | cut -d';' -f1
}


### SUBMIT ###

echo "Submitting pipeline (mode: $MODE, dependency type: $DEP_TYPE)"
echo "-------------------------------------------------------------"

SIM_DEP=""

# Main chain: network generation -> NCRG -> collection
if [[ "$MODE" == "all" || "$MODE" == "main" ]]; then
  J01=$(submit_sh 01_NetworkGeneratorComparison.sh)
  echo "01  Network generator comparison : $J01"

  J01a=$(submit_sh 01a_Run_NCRG.sh "$J01")
  echo "01a NCRG                         : $J01a  (after $J01)"

  J02=$(submit_r 02_collectGraphs.R collectGraphs "$J01a")
  echo "02  Collect graphs               : $J02  (after $J01a)"

  SIM_DEP="$J02"
fi

# Epidemic simulations: all four run in parallel once graphs are collected
if [[ "$MODE" == "all" || "$MODE" == "main" || "$MODE" == "sims" ]]; then
  J03a=$(submit_sh 03a_SEIR_Randomized.sh "$SIM_DEP")
  J03b=$(submit_sh 03b_SEIR_Factorial.sh  "$SIM_DEP")
  J04a=$(submit_sh 04a_SIS_Randomized.sh  "$SIM_DEP")
  J04b=$(submit_sh 04b_SIS_Factorial.sh   "$SIM_DEP")
  echo "03a SEIR randomized              : $J03a${SIM_DEP:+  (after $SIM_DEP)}"
  echo "03b SEIR factorial               : $J03b${SIM_DEP:+  (after $SIM_DEP)}"
  echo "04a SIS randomized               : $J04a${SIM_DEP:+  (after $SIM_DEP)}"
  echo "04b SIS factorial                : $J04b${SIM_DEP:+  (after $SIM_DEP)}"
fi

# Additional locations: independent of the main chain
if [[ "$MODE" == "all" || "$MODE" == "other" ]]; then
  J05=$(submit_sh 05_Graphs_Other_Locations.sh)
  echo "05  Graphs other locations       : $J05"

  J06=$(submit_r 06_collectGraphsOtherLocations.R collectOtherLoc "$J05")
  echo "06  Collect other locations      : $J06  (after $J05)"
fi

echo "-------------------------------------------------------------"
echo "Monitor with:  squeue -u \$USER"
echo "Cancel all:    scancel -u \$USER"
