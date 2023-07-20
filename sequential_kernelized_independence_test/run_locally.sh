#!/bin/bash
#SBATCH --mail-type=end,fail
#SBATCH --job-name="independence"
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=2000M
#SBATCH --array=1-3000

#### Your shell commands below this line ####

# simulate data

export SLURM_ARRAY_TASK_ID='1'
R CMD BATCH --no-save --no-restore simulate_data_skit.R

# run tests
export lmbd_type='ONS'
python3 simulations_skit.py
export lmbd_type='aGRAPA'
python3 simulations_skit.py

# remove data
rm "simulation_examples_$SLURM_ARRAY_TASK_ID.csv"
