#!/bin/bash

### MODULES

ml purge
ml R/4.2.1-foss-2022a # Adapt to your cluster's R version

### PARAMETERS

# Name of corresponding R script > don't change
R_SCRIPT="03a_SEIR_Randomized.R"  

# Master seed for reproducibility of manuscript results > don't change
MASTER_SEED=100

# Total number of tasks > you can change this based on how many simulations you want to run
# The current 10,000 tasks combined with 300 simulations per task yields a total of
# 3,000,000 simulations, resulting in 500,000 simulations per graph type (empirical graph incl.)
# TOTAL_TASKS and SIMS_PER_TASK thus function together.
# This should also correspond to #SBATCH --array
TOTAL_TASKS=10000

# How many simulations per task
SIMS_PER_TASK=300 # this is the current maximum given beta and number of seeds

# Infectious period > don't change for manuscript re-run
# Note that you could turn this into a range, but the script currently expects only a single value
# If you want to run a sensitivity analysis, you have to adapt the .R script too. 
DELTA=2

# Latent period > don't change for manuscript re-run
# Same comment as above
SIGMA=30

# Number of graphs per graph type > we selected 100 out of 5000 graphs. 
# Only change this if you changed that number
N_GRAPHS=100

# This variable corresponds to CURR_TYPE in the .R script
# In its current set-up, it searches for the Chad graphs in "Seed_5000" and
# each separate location too. You can change this by removing or adding
# NOTE: this needs to perfectly correspond with lines 406-416 in the .R file
LOCATIONS=("Seed_5000" "Romana" "Sabaneta" "Habi" "Hepang")

### RUN

# We run the total number of arrays (currently set at 10,000) for each location via a loop
for LOC in "${LOCATIONS[@]}"
do
  echo "Submitting 10,000 array tasks for location: $LOC"
  
  # Only here do we add the sbatch information 
  # Adapt CPU and memory as needed
  # Each separate job is currently put in a 30 minute queue (they are usually less busy) and
  # runs for 30 minutes maximum
  sbatch --job-name="SEIR_${LOC}" \
         --time=30:00 \
         --qos=30min \
         --mem-per-cpu=16G \
         --cpus-per-task=4 \
         --array=0-9999%1000 \
         --export=ALL,CURR_LOCATION="$LOC" \
         --wrap="Rscript $R_SCRIPT $MASTER_SEED $TOTAL_TASKS $SIMS_PER_TASK $DELTA $SIGMA $N_GRAPHS '$LOC'"
done

# Check if task completed or not
if [ $? -ne 0 ]; then
  echo "Task $SLURM_ARRAY_TASK_ID FAILED"
  exit 1
fi

echo "Task $SLURM_ARRAY_TASK_ID completed successfully"

### End of script