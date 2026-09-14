#!/bin/bash
#SBATCH --job-name=NetGen_New
#SBATCH --output="/scicore/home/chitnis/derkx0000/GraphComparison/SLURM/%x.%A_%a.out"
#SBATCH --error="/scicore/home/chitnis/derkx0000/GraphComparison/SLURM/%x.%A_%a.err"
#SBATCH --array=1-5000             # This creates 100 tasks
#SBATCH --cpus-per-task=1          # Each task is single-threaded
#SBATCH --mem=8G                   # Adjust based on your N_nodes
#SBATCH --time=06:00:00            # Adjust based on grid_optimization size

# Load R module (name depends on your cluster, e.g., R/4.2.0)
ml purge
ml R/4.2.1-foss-2022a # Load R version 4.2.1

# Define your combinations
COUNTRIES=("Guatemala" "Guatemala" "Indonesia" "Indonesia")
LOCATIONS=("Romana" "Sabaneta" "Habi" "Hepang" )
FILE_DIR="/scicore/home/chitnis/derkx0000/BaseData"
R_SCRIPT="06_Graphs_Other_Locations.R"

# The Seed is the Array Task ID
SEED=$SLURM_ARRAY_TASK_ID
echo "Processing Seed: $SEED"

for i in {0..3}
do
    CURRENT_COUNTRY=${COUNTRIES[$i]}
    CURRENT_LOCATION=${LOCATIONS[$i]}
    
    echo "Running for: $CURRENT_COUNTRY - $CURRENT_LOCATION (Seed: $SEED)"
    
    # Run the R script
    # IMPORTANT: Ensure 9_Graphs_Other_Locations.R saves its output with BOTH 
    # the Seed and the Location in the filename to prevent overwriting.
    Rscript $R_SCRIPT "$SEED" "$CURRENT_COUNTRY" "$CURRENT_LOCATION" "$FILE_DIR"
done