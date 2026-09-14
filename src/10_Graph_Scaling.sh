#!/bin/bash
#SBATCH --job-name=NetScal
#SBATCH --output="/scicore/home/chitnis/derkx0000/GraphComparison/SLURM/%x.%A_%a.out"
#SBATCH --error="/scicore/home/chitnis/derkx0000/GraphComparison/SLURM/%x.%A_%a.err"
#SBATCH --array=0-4999             # This creates 100 tasks
#SBATCH --cpus-per-task=1          # Each task is single-threaded
#SBATCH --mem=8G                   # Adjust based on your N_nodes
#SBATCH --time=06:00:00            # Adjust based on grid_optimization size

# Load R module (name depends on your cluster, e.g., R/4.2.0)
ml purge
ml R/4.2.1-foss-2022a # Load R version 4.2.1

# Run 10000 runs simultaneously
OFFSET=${OFFSET:-0}
ACTUAL_ID=$((SLURM_ARRAY_TASK_ID + OFFSET))
R_SCRIPT="10_Graph_Scaling.R"
SQUARE_EDGE=4
K_MAX=20
MAX_SEED=5000
echo "Processing Task ID: $ACTUAL_ID"

# Run R passing the Array Task ID as an argument
Rscript $R_SCRIPT $ACTUAL_ID $SQUARE_EDGE $K_MAX $MAX_SEED