#!/bin/bash
#SBATCH --job-name=00_NetworkGeneration
#SBATCH --output="/path/%x.%A_%a.out"
#SBATCH --error="/path/%x.%A_%a.err"
#SBATCH --array=1-5000%100         # This creates 100 simultaneous tasks in an array of 5000
#SBATCH --cpus-per-task=1          # Each task is single-threaded
#SBATCH --mem=8G                   # Adjust based on your N_nodes
#SBATCH --time=06:00:00            # Adjust based on grid_optimization size

### MODULES

ml purge
ml R/4.2.1-foss-2022a # Adapt to your cluster's R version

# Set task id and script
OFFSET=${OFFSET:-0}
ACTUAL_ID=$((SLURM_ARRAY_TASK_ID + OFFSET))
R_SCRIPT="01_NetworkGeneratorComparison.R"

### RUN

echo "Processing Task ID: $ACTUAL_ID"

# Run R passing the Array Task ID as an argument
Rscript $R_SCRIPT $ACTUAL_ID