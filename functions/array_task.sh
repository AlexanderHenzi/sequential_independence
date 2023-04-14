#!/bin/bash
#SBATCH --mail-type=end,fail
#SBATCH --job-name="independence"
#SBATCH --time=0:20:00
#SBATCH --mem-per-cpu=2000M
#SBATCH --partition=epyc2
#SBATCH --array=1-5000

#### Your shell commands below this line ####

module load R
R CMD BATCH --no-save --no-restore simulations_cluster.R

