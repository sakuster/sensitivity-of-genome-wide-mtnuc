#!/bin/bash

#SBATCH --time=36:00:00
#SBATACH --cpus-per-task=24
#SBATCH --mail-type=end
#SBATCH --mail-user=shady.kuster@colostate.edu
#SBATCH --output=nonnmt_HCXcorv4_10-20-2025

#activate R
source /home/skuster/anaconda3/etc/profile.d/conda.sh
conda activate rforme

#run R script
Rscript HUEX_STAC_non_nmt_allelefreq_calculationv4.R

#deactivate envs
conda deactivate

