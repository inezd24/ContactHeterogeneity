#!/bin/bash
#SBATCH --job-name=Newman_Patch
#SBATCH --output="/scicore/home/chitnis/derkx0000/GraphComparison/SLURM/%x.%A_%a.out"
#SBATCH --error="/scicore/home/chitnis/derkx0000/GraphComparison/SLURM/%x.%A_%a.err"
#SBATCH --array=1-5000
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00

# Load R module (name depends on your cluster, e.g., R/4.2.0)
ml purge
ml R/4.2.1-foss-2022a # Load R version 4.2.1

# Script requirements
OFFSET=${OFFSET:-0}
ACTUAL_ID=$((SLURM_ARRAY_TASK_ID + OFFSET))
R_SCRIPT="01a_Run_NCRG.R"

# Run script
echo "Processing Task ID: $ACTUAL_ID"
Rscript $R_SCRIPT $ACTUAL_ID