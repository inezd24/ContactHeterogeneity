#!/bin/bash

#SBATCH --job-name=simulation_plotting
#SBATCH --cpus-per-task=2 # change to 1 > only 1+ if parallelization
#SBATCH --mem=256G   # you will get entire node 
#SBATCH --time=2:00:00
#SBATCH --output="/scicore/home/chitnis/derkx0000/GraphComparison/SLURM/%x.%A_%a.out"
#SBATCH --error="/scicore/home/chitnis/derkx0000/GraphComparison/SLURM/%x.%A_%a.err" 
#SBATCH --array=0-3

# Define your combinations
MODELS=("SEIR" "SIS")
TYPE=("Factorial" "Factorial")
MASTER=100
SIMS=("23400000" "163800000")

# Map the array ID to the specific model/seed combo
MODEL=${MODELS[$SLURM_ARRAY_TASK_ID]}
TYPE=${TYPE[$SLURM_ARRAY_TASK_ID]}
SIM=${SIMS[$SLURM_ARRAY_TASK_ID]}
RSCRIPT="05_Graph_Plotting.R"

# Load R module (adjust name based on your cluster)
ml purge
ml R/4.2.1-foss-2022a # Load R version 4.2.1

echo "Processing Combo: $MODEL model with $MASTER seeds with $SIM total simulations for the $TYPE location."

# Execute the R script
Rscript $RSCRIPT $MODEL $TYPE $MASTER $SIM