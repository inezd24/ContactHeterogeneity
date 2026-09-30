#!/bin/bash
#SBATCH --job-name=Net_OtherLoc
#SBATCH --output="/path/%x.%A_%a.out"
#SBATCH --error="/path/%x.%A_%a.err"
#SBATCH --array=1-5000            
#SBATCH --cpus-per-task=1  
#SBATCH --mem=8G             
#SBATCH --time=06:00:00       

### MODULES

ml purge
ml R/4.2.1-foss-2022a # Adapt to your cluster's R version

### PARAMETERS

# Which countries we're looking at (for the data files)
# The repetitions need to match the locations below
COUNTRIES=("Guatemala" "Guatemala" "Indonesia" "Indonesia")

# Which locations (needs to match countries)
LOCATIONS=("Romana" "Sabaneta" "Habi" "Hepang" )

# Directory of base folder where the country_location folders are
FILE_DIR="/scicore/home/chitnis/derkx0000/BaseData"

# Corresponding R script
R_SCRIPT="06_Graphs_Other_Locations.R"

# The Seed is the Array Task ID
SEED=$SLURM_ARRAY_TASK_ID

### RUN

echo "Processing Seed: $SEED"

# We loop over the four location/country combinations
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