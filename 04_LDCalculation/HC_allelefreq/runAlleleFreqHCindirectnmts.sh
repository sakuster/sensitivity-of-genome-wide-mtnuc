#!/bin/bash

#SBATCH --time=2:00:00
#SBATACH --cpus-per-task=24
#SBATCH --mail-type=end
#SBATCH --mail-user=shady.kuster@colostate.edu
#SBATCH --output=indirectnmt_HC_Xcor_2025-03-18

#activate R
source /home/skuster/anaconda3/etc/profile.d/conda.sh
conda activate rforme

#run R script
Rscript HUEX_STAC_indirectn-mt_allelefreq_calculation.R

#deactivate envs
conda deactivate

