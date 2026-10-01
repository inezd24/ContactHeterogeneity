#!/bin/bash
#SBATCH --job-name=SEIR_fact
#SBATCH --time=30:00
#SBATCH --qos=30min
#SBATCH --mem-per-cpu=16G
#SBATCH --cpus-per-task=4
#SBATCH --array=0-9999
#SBATCH --output="/path/%x.%A_%a.out"
#SBATCH --error="/path/%x.%A_%a.err"

### MODULES

ml purge
ml R/4.2.1-foss-2022a # Adapt to your cluster's R version

### PARAMETERS

# Change according to where your data is stored > change for your own path
GRAPH_FILE="/scicore/home/chitnis/derkx0000/Derkx2026_publication/Net_Sens/Graphs/Seed_5000/MASTER_ENSEMBLE_FOR_DM.rds"

# Name of corresponding R script > don't change
R_SCRIPT="03b_SEIR_Factorial.R"

# Master seed for reproducibility of manuscript results > don't change
MASTER_SEED=100

# Total number of tasks > you can change this based on how many simulations you want to run
# This should correspond to #SBATCH --array. For manuscript, replication, use 10,000
TOTAL_TASKS=10000

# Number of graphs per graph type > we selected 100 out of 5000 graphs. 
# Only change this if you changed that number
N_GRAPHS=100

# Infectious period > don't change for manuscript re-run
# Note that you could turn this into a range, but the script currently expects only a single value
# If you want to run a sensitivity analysis, you have to adapt the .R script too. 
DELTA=2

# Latent period > don't change for manuscript re-run
# Same comment as above
SIGMA=30

# Number of repetitions for the empirical graph
N_REPS=1000


### RUN

echo "Starting task $SLURM_ARRAY_TASK_ID"

Rscript $R_SCRIPT \
    "$GRAPH_FILE" \
    $MASTER_SEED \
    $TOTAL_TASKS \
    $DELTA \
    $SIGMA \
    $N_GRAPHS \
    $N_REPS

### CHECK

if [ $? -ne 0 ]; then
  echo "Task $SLURM_ARRAY_TASK_ID FAILED"
  exit 1
fi

echo "Task $SLURM_ARRAY_TASK_ID completed successfully"