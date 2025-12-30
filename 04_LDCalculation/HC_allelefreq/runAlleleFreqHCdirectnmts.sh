#!/bin/bash

#SBATCH --time=00:10:00
#SBATACH --cpus-per-task=24
#SBATCH --mail-type=end
#SBATCH --mail-user=shady.kuster@colostate.edu
#SBATCH --output=v4_directnmt_mt1_Xcor_10-20-2025

#activate R
source /home/skuster/anaconda3/etc/profile.d/conda.sh
conda activate rforme

#run R script
Rscript HUEX_STAC_directn-mt_allelefreqv4_calculation.R

#deactivate envs
conda deactivate

