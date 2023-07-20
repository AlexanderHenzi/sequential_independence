#!/bin/bash
#SBATCH --mail-type=end,fail
#SBATCH --job-name="independence"
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=2000M
#SBATCH --array=1-3000

#### Your shell commands below this line ####

# modules
module load gcc/8.2.0
module load r/4.1.3
module load python/3.8.5

# simulate data
R CMD BATCH --no-save --no-restore simulate_data_skit.R

# run tests
export lmbd_type='ONS'
python3 simulations_skit.py
export lmbd_type='aGRAPA'
python3 simulations_skit.py

# remove data
rm "simulation_examples_$SLURM_ARRAY_TASK_ID.csv"
