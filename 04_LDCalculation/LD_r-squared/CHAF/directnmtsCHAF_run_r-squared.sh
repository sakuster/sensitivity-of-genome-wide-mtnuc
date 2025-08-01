#!/bin/bash

#SBATCH --time=00:10:00
#SBATACH --cpus-per-task=24
#SBATCH --mail-type=end
#SBATCH --mail-user=shady.kuster@colostate.edu
#SBATCH --output=directnmtCHAF_rsquared_mt1_08-29-24

#activate R
source /home/skuster/anaconda3/etc/profile.d/conda.sh
conda activate rforme

#run R script
Rscript directnmt_rsquared_calculation_CHAF.R

#deactivate envs
conda deactivate

