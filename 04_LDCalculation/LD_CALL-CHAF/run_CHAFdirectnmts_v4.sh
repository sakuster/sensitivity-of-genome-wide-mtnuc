#!/bin/bash

#SBATCH --time=00:10:00
#SBATCH --mail-type=end
#SBATCH --mail-user=shady.kuster@colostate.edu
#SBATCH --output=v4_directnmt_mt1_CHAF_10-16-2025

#activate R
source /home/skuster/anaconda3/etc/profile.d/conda.sh
conda activate rforme

#run R script
Rscript directnmt_LD_calculation_v4_CHAF.R

#deactivate envs
conda deactivate

