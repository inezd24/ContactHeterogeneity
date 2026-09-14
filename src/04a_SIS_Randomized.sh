#!/bin/bash
#SBATCH --job-name=SIS_rand
#SBATCH --time=30:00
#SBATCH --qos=30min
#SBATCH --mem-per-cpu=32G
#SBATCH --cpus-per-task=4
#SBATCH --array=0-999
#SBATCH --output="/scicore/home/chitnis/derkx0000/GraphComparison/SLURM/%x.%A_%a.out"
#SBATCH --error="/scicore/home/chitnis/derkx0000/GraphComparison/SLURM/%x.%A_%a.err"

### MODULES
ml purge
ml R/4.2.1-foss-2022a

### PARAMETERS
GRAPH_FILE="/scicore/home/chitnis/derkx0000/GraphComparison/Net_Sens/Graphs/Seed_5000/MASTER_ENSEMBLE_FOR_DM.rds"
R_SCRIPT="04a_SIS_Randomized.R"
MASTER_SEED=100
TOTAL_TASKS=1000
SIMS_PER_TASK=6000
N_GRAPHS=100

### RUN R SCRIPT
echo "Starting task $SLURM_ARRAY_TASK_ID"
Rscript $R_SCRIPT \
  "$GRAPH_FILE" \
  $MASTER_SEED \
  $TOTAL_TASKS \
  $SIMS_PER_TASK \
  $N_GRAPHS

if [ $? -ne 0 ]; then
  echo "Task $SLURM_ARRAY_TASK_ID FAILED"
  exit 1
fi

echo "Task $SLURM_ARRAY_TASK_ID completed successfully"
