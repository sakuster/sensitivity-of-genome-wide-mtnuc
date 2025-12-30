#!/bin/bash

#SBATCH --time=00:10:00
#SBATACH --cpus-per-task=24
#SBATCH --mail-type=end
#SBATCH --mail-user=shady.kuster@colostate.edu
#SBATCH --output=v4_directnmt_mt1_CALL_10-15-2025

#activate R
source /home/skuster/anaconda3/etc/profile.d/conda.sh
conda activate rforme

#run R script
Rscript directnmt_LD_calculation_v4_CALL.R

#deactivate envs
conda deactivate

