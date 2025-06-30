#!/bin/bash

#SBATCH --time=36:00:00
#SBATACH --cpus-per-task=24
#SBATCH --mail-type=end
#SBATCH --mail-user=shady.kuster@colostate.edu
#SBATCH --output=v3_nonnmt_CALLmt1_08-28-24

#activate R
source /home/skuster/anaconda3/etc/profile.d/conda.sh
conda activate rforme

#run R script
Rscript LD_calculation_v3_CALLmt1nonnmt.R

#deactivate envs
conda deactivate

