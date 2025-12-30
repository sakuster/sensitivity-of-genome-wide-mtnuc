#!/bin/bash

#SBATCH --time=2:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=shady.kuster@colostate.edu
#SBATCH --output=v4_indirectnmt_CHAF_10-16-2025

#activate R
source /home/skuster/anaconda3/etc/profile.d/conda.sh
conda activate rforme

#run R script
Rscript LD_calculation_v4_CHAFmt1indir.R

#deactivate envs
conda deactivate

