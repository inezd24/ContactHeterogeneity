#!/bin/bash
#SBATCH --job-name=Newman_Patch
#SBATCH --output="/path/%x.%A_%a.out"
#SBATCH --error="/path/%x.%A_%a.err"
#SBATCH --array=1-5000
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00

### MODULES

ml purge
ml R/4.2.1-foss-2022a # Adapt to your cluster's R version

### PARAMETERS

# Offset for task id (you can ignore)
OFFSET=${OFFSET:-0}

# Task id
ACTUAL_ID=$((SLURM_ARRAY_TASK_ID + OFFSET))

# Script name
R_SCRIPT="01a_Run_NCRG.R"

# Run script
echo "Processing Task ID: $ACTUAL_ID"
Rscript $R_SCRIPT $ACTUAL_ID