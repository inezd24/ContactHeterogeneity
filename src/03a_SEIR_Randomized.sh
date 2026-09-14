#!/bin/bash

### MODULES
ml purge
ml R/4.2.1-foss-2022a

### PARAMETERS
R_SCRIPT="03a_SEIR_Randomized.R"
MASTER_SEED=100
TOTAL_TASKS=10000
SIMS_PER_TASK=300 # this is the current maximum given beta and number of seeds
DELTA=2
SIGMA=30
N_GRAPHS=100

# Define your 4 target locations
LOCATIONS=("Romana" "Sabaneta" "Habi" "Hepang")

# Loop through and launch a separate 10,000-job array for each location
for LOC in "${LOCATIONS[@]}"
do
  echo "Submitting 10,000 array tasks for location: $LOC"
  
  sbatch --job-name="SEIR_${LOC}" \
         --time=30:00 \
         --qos=30min \
         --mem-per-cpu=16G \
         --cpus-per-task=4 \
         --array=0-9999%1000 \
         --export=ALL,CURR_LOCATION="$LOC" \
         --wrap="Rscript $R_SCRIPT $MASTER_SEED $TOTAL_TASKS $SIMS_PER_TASK $DELTA $SIGMA $N_GRAPHS '$LOC'"
done

if [ $? -ne 0 ]; then
  echo "Task $SLURM_ARRAY_TASK_ID FAILED"
  exit 1
fi

echo "Task $SLURM_ARRAY_TASK_ID completed successfully"