#!/bin/bash
#SBATCH --job-name=ProcGraphs
#SBATCH --output="/scicore/home/chitnis/derkx0000/GraphComparison/SLURM/%x.%A_%a.out"
#SBATCH --error="/scicore/home/chitnis/derkx0000/GraphComparison/SLURM/%x.%A_%a.err"
#SBATCH --cpus-per-task=4          # Each task is single-threaded
#SBATCH --mem=128G                   # Adjust based on your N_nodes
#SBATCH --time=00:30:00            # Adjust based on grid_optimization size

# Load R module (name depends on your cluster, e.g., R/4.2.0)
ml purge
ml R/4.2.1-foss-2022a # Load R version 4.2.1

# Script
R_SCRIPT="11_Process_Scaled_Graphs.R"

# Run R 
Rscript $R_SCRIPT 